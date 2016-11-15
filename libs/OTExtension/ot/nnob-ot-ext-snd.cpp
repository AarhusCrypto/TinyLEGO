/*
 * nnob-ot-ext-snd.cpp
 *
 *  Created on: Mar 23, 2015
 *      Author: mzohner
 */


#include "nnob-ot-ext-snd.h"

BOOL NNOBOTExtSnd::sender_routine(uint32_t id, uint64_t myNumOTs) {
	uint64_t myStartPos = id * myNumOTs;
	uint64_t wd_size_bits = m_nBlockSizeBits;
	uint64_t processedOTBlocks = min((uint64_t) NUMOTBLOCKS, ceil_divide(myNumOTs, wd_size_bits));
	uint64_t OTsPerIteration = processedOTBlocks * wd_size_bits;
	uint64_t tmpctr, tmpotlen;
	uint32_t nchans = 2;
	bool use_mat_chan = (m_eSndOTFlav == Snd_GC_OT || m_bUseMinEntCorRob);
	if(use_mat_chan) {
		nchans = 3;
	}

	channel* ot_chan = new channel(OT_BASE_CHANNEL+nchans*id, m_cRcvThread, m_cSndThread);
	channel* check_chan = new channel(OT_BASE_CHANNEL+nchans*id + 1, m_cRcvThread, m_cSndThread);
	channel* mat_chan;
	if(use_mat_chan) {
		mat_chan = new channel(nchans*id+2, m_cRcvThread, m_cSndThread);
	}

	myNumOTs = min(myNumOTs + myStartPos, m_nOTs) - myStartPos;
	uint64_t lim = myStartPos + myNumOTs;
	uint64_t** rndmat;

	// The vector with the received bits
	CBitVector vRcv(m_nBaseOTs * OTsPerIteration);
	vRcv.Reset();

	// Holds the reply that is sent back to the receiver
	uint32_t numsndvals = 2;
	CBitVector* vSnd;

	CBitVector* seedbuf = new CBitVector[m_nSndVals];
	for (uint32_t u = 0; u < m_nSndVals; u++)
		seedbuf[u].Create(OTsPerIteration * m_cCrypt->get_aes_key_bytes() * 8);
#ifdef ZDEBUG
	cout << "seedbuf size = " <<OTsPerIteration * AES_KEY_BITS << endl;
#endif
	vSnd = new CBitVector[numsndvals];
	for (uint32_t i = 0; i < numsndvals; i++) {
		vSnd[i].Create(OTsPerIteration * m_nBitLength);
	}

	// Contains the parts of the V matrix
	CBitVector Q(wd_size_bits * OTsPerIteration);
	mask_buf_t tmpmaskbuf;

	uint64_t OT_ptr = myStartPos;

	uint8_t *rcvbuftmpptr, *rcvbufptr;

	queue<nnob_snd_check_t*> check_queue;
	queue<mask_buf_t> mask_queue;

	uint32_t startpos = 0;
	if(m_eRecOTFlav==Rec_R_OT) {
		startpos = 1;
	}

	if(m_eSndOTFlav == Snd_GC_OT) {
		initRndMatrix(&rndmat, m_nBitLength, m_nBaseOTs);
	}

#ifdef OTTiming
	double totalMtxTime = 0, totalTnsTime = 0, totalHshTime = 0, totalRcvTime = 0, totalSndTime = 0, totalUnMaskTime = 0,
			totalHashCheckTime = 0;
	timeval tempStart, tempEnd;
#endif

	while (OT_ptr < lim) //do while there are still transfers missing
	{
		processedOTBlocks = min((uint64_t) NUMOTBLOCKS, ceil_divide(lim - OT_ptr, wd_size_bits));
		OTsPerIteration = processedOTBlocks * wd_size_bits;

#ifdef ZDEBUG
		cout << "Processing block " << nProgress << " with length: " << OTsPerIteration << ", and limit: " << lim << endl;
#endif

#ifdef OTTiming
		gettimeofday(&tempStart, NULL);
#endif
		rcvbufptr = ot_chan->blocking_receive_id_len(&rcvbuftmpptr, &tmpctr, &tmpotlen);
		//vRcv.AttachBuf(rcvbuftmpptr, bits_in_bytes(m_nBaseOTs * OTsPerIteration));
		vRcv.SetBytes(rcvbuftmpptr, bits_in_bytes(OTsPerIteration*startpos), bits_in_bytes((m_nBaseOTs-startpos)*OTsPerIteration));
		free(rcvbufptr);
		//vRcv.PrintHex();
#ifdef OTTiming
		gettimeofday(&tempEnd, NULL);
		totalRcvTime += getMillies(tempStart, tempEnd);
		gettimeofday(&tempStart, NULL);
#endif
		BuildQMatrix(&Q, OT_ptr, processedOTBlocks, m_tBaseOTKeys.front());
#ifdef OTTiming
		gettimeofday(&tempEnd, NULL);
		totalMtxTime += getMillies(tempStart, tempEnd);
		gettimeofday(&tempStart, NULL);
#endif
		check_queue.push(UpdateCheckBuf(Q.GetArr(), vRcv.GetArr(), OT_ptr, processedOTBlocks, check_chan));
		FillAndSendRandomMatrix(rndmat, mat_chan);
#ifdef OTTiming
		gettimeofday(&tempEnd, NULL);
		totalHashCheckTime += getMillies(tempStart, tempEnd);
		gettimeofday(&tempStart, NULL);
#endif
		UnMaskBaseOTs(&Q, &vRcv, m_tBaseOTChoices.front(), processedOTBlocks);

		GenerateSendAndXORCorRobVector(&Q, OTsPerIteration, mat_chan);

#ifdef OTTiming
		gettimeofday(&tempEnd, NULL);
		totalUnMaskTime += getMillies(tempStart, tempEnd);
		gettimeofday(&tempStart, NULL);
#endif
		Q.Transpose(wd_size_bits, OTsPerIteration);
#ifdef OTTiming
		gettimeofday(&tempEnd, NULL);
		totalTnsTime += getMillies(tempStart, tempEnd);
		gettimeofday(&tempStart, NULL);
#endif
		HashValues(&Q, seedbuf, vSnd, m_tBaseOTChoices.front(), OT_ptr, min(lim - OT_ptr, OTsPerIteration), rndmat);
#ifdef OTTiming
		gettimeofday(&tempEnd, NULL);
		totalHshTime += getMillies(tempStart, tempEnd);
		gettimeofday(&tempStart, NULL);
#endif

		//TODO: outsource into method
		tmpmaskbuf.otid = OT_ptr;
		tmpmaskbuf.otlen = min(lim - OT_ptr, OTsPerIteration);
		tmpmaskbuf.maskbuf = new CBitVector[numsndvals];
		for(uint32_t i = 0; i < numsndvals; i++)
			tmpmaskbuf.maskbuf[i].Copy(vSnd[i]);
		mask_queue.push(tmpmaskbuf);

		if(check_chan->data_available()) {
			assert(CheckConsistency(&check_queue, check_chan));
			tmpmaskbuf = mask_queue.front();
			mask_queue.pop();
			MaskAndSend(tmpmaskbuf.maskbuf, tmpmaskbuf.otid, tmpmaskbuf.otlen, ot_chan);
			for(uint32_t i = 0; i < numsndvals; i++)
				tmpmaskbuf.maskbuf[i].delCBitVector();
			delete tmpmaskbuf.maskbuf;
		}
#ifdef OTTiming
		gettimeofday(&tempEnd, NULL);
		totalSndTime += getMillies(tempStart, tempEnd);
#endif
		OT_ptr += min(lim - OT_ptr, OTsPerIteration);
		Q.Reset();

		//free(rcvbufptr);
	}

	while(!check_queue.empty()) {
		if(check_chan->data_available()) {
			assert(CheckConsistency(&check_queue, check_chan));//TODO assert
			tmpmaskbuf = mask_queue.front();
			mask_queue.pop();
			MaskAndSend(tmpmaskbuf.maskbuf, tmpmaskbuf.otid, tmpmaskbuf.otlen, ot_chan);
			for(uint32_t i = 0; i < numsndvals; i++)
				tmpmaskbuf.maskbuf[i].delCBitVector();
			delete tmpmaskbuf.maskbuf;
		}
	}

	ot_chan->synchronize_end();
	check_chan->synchronize_end();

	Q.delCBitVector();
	vRcv.delCBitVector();

	for (uint32_t u = 0; u < m_nSndVals; u++)
		seedbuf[u].delCBitVector();
#ifndef ABY_OT
	delete seedbuf;
#endif
	for (uint32_t i = 0; i < numsndvals; i++)
		vSnd[i].delCBitVector();
#ifndef ABY_OT
	if (numsndvals > 0)
		delete vSnd;
#endif

	if(use_mat_chan) {
		mat_chan->synchronize_end();
	}

	if(m_eSndOTFlav == Snd_GC_OT) {
		freeRndMatrix(rndmat, m_nBaseOTs);
	}
#ifdef OTTiming
	cout << "Sender time benchmark for performing " << myNumOTs << " OTs on " << m_nBitLength << " bit strings" << endl;
	cout << "Time needed for: " << endl;
	cout << "\t Matrix Generation:\t" << totalMtxTime << " ms" << endl;
	cout << "\t BaseOT Unmasking:\t" << totalUnMaskTime << " ms" << endl;
	cout << "\t Check Hashing:\t" << totalHashCheckTime << " ms" << endl;
	cout << "\t Sending Matrix:\t" << totalSndTime << " ms" << endl;
	cout << "\t Transposing Matrix:\t" << totalTnsTime << " ms" << endl;
	cout << "\t Hashing Matrix:\t" << totalHshTime << " ms" << endl;
	cout << "\t Receiving Values:\t" << totalRcvTime << " ms" << endl;
#endif


	delete ot_chan;
	delete check_chan;
	return TRUE;
}


void NNOBOTExtSnd::FillAndSendRandomMatrix(uint64_t **rndmat, channel* mat_chan) {
	if(m_eSndOTFlav == Snd_GC_OT) {
		uint8_t* rnd_seed = (uint8_t*) malloc(m_nSymSecParam);
		m_cCrypt->gen_rnd(rnd_seed, m_nSymSecParam);
		mat_chan->send(rnd_seed, m_nSymSecParam);
		fillRndMatrix(rnd_seed, rndmat, m_nBitLength, m_nBaseOTs, m_cCrypt);
		free(rnd_seed);
	}
}



nnob_snd_check_t* NNOBOTExtSnd::UpdateCheckBuf(uint8_t* tocheckseed, uint8_t* tocheckrcv, uint64_t otid, uint64_t numblocks, channel* check_chan) {
	uint64_t rowbytelen = m_nBlockSizeBytes * numblocks;
	uint64_t checkbytelen = min(rowbytelen, bits_in_bytes(m_nOTs - otid));

	uint8_t* hash_buf = (uint8_t*) malloc(SHA512_DIGEST_LENGTH);
	uint8_t* tmpbuf = (uint8_t*) malloc(rowbytelen);
	uint8_t *idaptr, *idbptr;
	nnob_snd_check_t* check_buf = (nnob_snd_check_t*) malloc(sizeof(nnob_snd_check_t));
	//check_buf.rcv_chk_buf = (uint8_t*) malloc(m_nChecks * OWF_BYTES);
	check_buf->chk_buf = (uint8_t*) malloc(m_nChecks * OWF_BYTES);

	uint8_t *chk_buf_ptr = check_buf->chk_buf;
	uint8_t *idatmpbuf = (BYTE*) malloc(sizeof(BYTE) * rowbytelen);
	uint8_t *idbtmpbuf = (BYTE*) malloc(sizeof(BYTE) * rowbytelen);
	uint8_t *seedptr, *rcvptr;

	//uint32_t blockoffset = ceil_divide(otid, NUMOTBLOCKS * m_nBlockSizeBits);
	uint32_t blockid = 0; //TODO bring in as soon as 3-step OT is implemented

	check_buf->otid = otid;
	check_buf->numblocks = numblocks;
	check_buf->perm = (linking_t*) malloc(sizeof(linking_t*) * m_nChecks);
	check_buf->permchoicebits = (BYTE*) malloc(sizeof(uint8_t) * m_nChecks);

	genRandomMapping(check_buf->perm, m_nBaseOTs);

	for(uint32_t i = 0; i < m_nChecks; i++) {
		check_buf->permchoicebits[i] = m_tBaseOTChoices.front()->GetBit(blockid * m_nBaseOTs + check_buf->perm[i].ida) ^
				m_tBaseOTChoices.front()->GetBit(blockid * m_nBaseOTs + check_buf->perm[i].idb);
	}

	//right now the checkbytelen needs to be a multiple of AES_BYTES
	assert(ceil_divide(rowbytelen, OWF_BYTES) * OWF_BYTES == rowbytelen);
#ifdef DEBUG_NNOB_CHECKS
	m_vU.PrintHex();
#endif

	for(uint64_t i = 0; i < m_nChecks; i++, chk_buf_ptr+=OWF_BYTES) {

		if(m_tBaseOTChoices.front()->GetBit(blockid * m_nBaseOTs + check_buf->perm[i].ida) == 0) {
			memcpy(idatmpbuf, tocheckseed + check_buf->perm[i].ida * rowbytelen, rowbytelen);
		} else {
			seedptr = tocheckseed + check_buf->perm[i].ida * rowbytelen;
			rcvptr = tocheckrcv + check_buf->perm[i].ida * rowbytelen;
			for(int j = 0; j < rowbytelen/sizeof(uint64_t); j++) {
				((uint64_t*) idatmpbuf)[j] = ((uint64_t*) seedptr)[j] ^ ((uint64_t*) rcvptr)[j];
			}
		}

		if(m_tBaseOTChoices.front()->GetBit(blockid * m_nBaseOTs + check_buf->perm[i].idb) == 0) {
			memcpy(idbtmpbuf, tocheckseed + check_buf->perm[i].idb * rowbytelen, rowbytelen);
		} else {
			seedptr = tocheckseed + check_buf->perm[i].idb * rowbytelen;
			rcvptr = tocheckrcv + check_buf->perm[i].idb * rowbytelen;
			for(int j = 0; j < rowbytelen/sizeof(uint64_t); j++) {
				((uint64_t*) idbtmpbuf)[j] = ((uint64_t*) seedptr)[j] ^ ((uint64_t*) rcvptr)[j];
			}
		}

	/*#ifdef DEBUG_NNOB_CHECKS
		cout << "seedA: " <<  (hex) << ((uint64_t*) (tocheckseed + check_buf.perm[i].ida * rowbytelen))[0] <<
				", rcvA: " << ((uint64_t*) (tocheckrcv + check_buf.perm[i].ida * rowbytelen))[0] << (dec) << endl;
		cout << "seedB: " <<  (hex) << ((uint64_t*) (tocheckseed + check_buf.perm[i].idb * rowbytelen))[0] <<
				", rcvB: " << ((uint64_t*) (tocheckrcv + check_buf.perm[i].idb * rowbytelen))[0] << (dec) << endl;
		cout << "input to owf " <<  (hex) << ((uint64_t*) idatmpbuf)[0] << ", " << ((uint64_t*) idbtmpbuf)[0] << (dec) << endl;
	#endif*/

		memset(tmpbuf, 0, rowbytelen);

		for(uint64_t j = 0; j < rowbytelen/sizeof(uint64_t); j++) {
			((uint64_t*) tmpbuf)[j] = ((uint64_t*) tmpbuf)[j] ^ ((uint64_t*) idatmpbuf)[j] ^ ((uint64_t*) idbtmpbuf)[j];
		}

#ifdef DEBUG_NNOB_CHECKS_INPUT
		cout << "XOR-OWF Input:\t" << (hex);
		for(uint32_t t = 0; t < checkbytelen; t++) {
			cout << setw(2) << setfill('0') << (uint32_t) tmpbuf[t];
		}
		cout << (dec) << endl;
#endif
#ifdef AES_OWF
			owf(&aesowfkey, rowbytelen, tmpbuf, resbuf);
#else
		//m_cCrypt->hash_buf(chk_buf_ptr, OWF_BYTES, tmpbuf, checkbytelen, hash_buf);//hash_buf, rowbytelen, tmpbuf, resbuf, hash_buf);
		sha512_hash(chk_buf_ptr, OWF_BYTES, tmpbuf, checkbytelen, hash_buf);
#endif
#ifdef DEBUG_NNOB_CHECKS_OUTPUT
		cout << "XOR-OWF Output:\t" << (hex);
		for(uint32_t t = 0; t < OWF_BYTES; t++) {
			cout << (uint32_t) chk_buf_ptr[t];
		}
		cout << (dec) << endl;
#endif
		//XORandOWF(idatmpbuf, idbtmpbuf,	checkbytelen, tmpbuf, chk_buf_ptr, hash_buf);
	}

/*	for(uint64_t i = 0; i < m_nChecks; i++, seedcheckbufptr+=OWF_BYTES, rcvcheckbufptr+=OWF_BYTES) {
		memset(tmpbuf, 0, rowbytelen);
#ifdef DEBUG_ALSZ_CHECKS
		cout << i << "-th check between " << check_buf.perm[i].ida << " and " << check_buf.perm[i].idb << ": " << endl;
#endif
		XORandOWF(tocheckseed + check_buf.perm[i].ida * rowbytelen, tocheckseed + check_buf.perm[i].idb * rowbytelen,
				rowbytelen, tmpbuf, seedcheckbufptr, hash_buf);
		XORandOWF(tocheckrcv + check_buf.perm[i].ida * rowbytelen, tocheckrcv + check_buf.perm[i].idb * rowbytelen,
				rowbytelen, tmpbuf, rcvcheckbufptr, hash_buf);
	}*/

	free(tmpbuf);
	free(hash_buf);
	free(idatmpbuf);
	free(idbtmpbuf);

	//Send the permutation and the XORed bits over to the receiver
	check_chan->send_id_len((uint8_t*) check_buf->perm, sizeof(linking_t) * m_nChecks, otid, numblocks);
	check_chan->send(check_buf->permchoicebits, m_nChecks);

	return check_buf;
}

void NNOBOTExtSnd::XORandOWF(uint8_t* idaptr, uint8_t* idbptr, uint64_t rowbytelen, uint8_t* tmpbuf,
		uint8_t* resbuf, uint8_t* hash_buf) {

	memset(tmpbuf, 0, rowbytelen);

	for(uint64_t j = 0; j < rowbytelen/sizeof(uint64_t); j++) {
		((uint64_t*) tmpbuf)[j] = ((uint64_t*) tmpbuf)[j] ^ ((uint64_t*) idaptr)[j] ^ ((uint64_t*) idbptr)[j];
	}

#ifdef DEBUG_NNOB_CHECKS_INPUT
	cout << "XOR-OWF Input:\t" << (hex);
	for(uint32_t t = 0; t < rowbytelen; t++) {
		cout << setw(2) << setfill('0') << (uint32_t) tmpbuf[t];
	}
	cout << (dec) << endl;
#endif
#ifdef AES_OWF
		owf(&aesowfkey, rowbytelen, tmpbuf, resbuf);
#else
	m_cCrypt->hash_buf(resbuf, OWF_BYTES, tmpbuf, rowbytelen, hash_buf);//hash_buf, rowbytelen, tmpbuf, resbuf, hash_buf);
#endif
#ifdef DEBUG_NNOB_CHECKS_OUTPUT
	cout << "XOR-OWF Output:\t" << (hex);
	for(uint32_t t = 0; t < OWF_BYTES; t++) {
		cout << (uint32_t) resbuf[t];
	}
	cout << (dec) << endl;
#endif
}

BOOL NNOBOTExtSnd::CheckConsistency(queue<nnob_snd_check_t*>* check_buf_q, channel* check_chan) {
	uint8_t *rcvhashbufptr, *rcvhashbuf;
	uint64_t tmpid, tmpnblocks;

	uint8_t* rcvbuf = check_chan->blocking_receive_id_len(&rcvhashbuf, &tmpid, &tmpnblocks);

	nnob_snd_check_t* check_buf = check_buf_q->front();
	check_buf_q->pop();

	//Should be fine since the blocks are handled sequentially - but recheck anyway
	assert(check_buf->otid == tmpid);
	assert(check_buf->numblocks == tmpnblocks);

	rcvhashbufptr = rcvhashbuf;

	//Very simple : just go over both arrays and check equality
	uint64_t* rcvbufptr = (uint64_t*) rcvhashbuf;
	uint64_t* chkbufptr = (uint64_t*) check_buf->chk_buf;
	for(uint32_t i = 0; i < m_nChecks; i++) {
		for(uint32_t j = 0; j < OWF_BYTES / sizeof(uint64_t); j++, rcvbufptr++, chkbufptr++) {
			if(*rcvbufptr != *chkbufptr) {
	#ifdef DEBUG_NNOB_CHECKS
				cout << "Error in " << i <<"-th consistency check: " << endl;
				cout << "Receiver hash = " << (hex) << *rcvbufptr << ", my hash: " << *chkbufptr << endl;
	#endif
				return FALSE;
			}
		}
	}

	//free the receive buffer
	free(rcvbuf);
	free(check_buf->perm);
	free(check_buf->chk_buf);
	free(check_buf->permchoicebits);
	free(check_buf);

	return TRUE;
}


void NNOBOTExtSnd::genRandomMapping(linking_t* outperm, uint32_t nids) {
	uint32_t nperms = nids / 2;
	uint32_t i;
	uint32_t tmpidx = 0;

	uint32_t* mapping = (uint32_t*) malloc(sizeof(uint32_t) * nperms);
	for(i = 0; i < nperms; i++) {
		mapping[i] = nperms + i;
	}

	//shuffle content randomly using Knuths permutation algorithm
	for(i = 0; i < nperms; i++) {
		m_cCrypt->gen_rnd_uniform(&tmpidx, nperms);
		swap(mapping[i], mapping[tmpidx]);
	}

	for(i = 0; i < nperms; i++) {
		outperm[i].ida = i;
		outperm[i].idb = mapping[i];
		assert(outperm[i].idb >= nperms && outperm[i].idb < nids);
		//cout << i << " checked against " << outperm[i].idb << endl;
	}

	free(mapping);
}


void NNOBOTExtSnd::ComputeBaseOTs(field_type ftype) {
	if(m_bDoBaseOTs) {
		m_cBaseOT = new SimpleOT(m_cCrypt, ftype);
		ComputePKBaseOTs();
		delete m_cBaseOT;
	} else {
		//recursive call
	}
}
