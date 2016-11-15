/*
 * ot-extension-receiver.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: mzohner
 */
#include "ot-ext-rec.h"

BOOL OTExtRec::receive(uint64_t numOTs, uint64_t bitlength, uint64_t nsndvals, CBitVector* choices, CBitVector* ret,
		snd_ot_flavor stype, rec_ot_flavor rtype, uint32_t numThreads, MaskingFunction* unmaskfct) {
	m_nOTs = numOTs;
	m_nBitLength = bitlength;
	m_nSndVals = nsndvals;
	m_vChoices = choices;
	m_vRet = ret;
	m_eSndOTFlav = stype;
	m_eRecOTFlav = rtype;
	m_fMaskFct = unmaskfct;

	assert(pad_to_power_of_two(m_nSndVals) == m_nSndVals);

	return start_receive(numThreads);
}
;

//Initialize and start numThreads OTSenderThread
BOOL OTExtRec::start_receive(uint32_t numThreads) {
	if (m_nOTs == 0)
		return true;

	if(numThreads * m_nBlockSizeBits > m_nOTs && numThreads > 1) {
		cerr << "Decreasing nthreads from " << numThreads << " to " << max(m_nOTs / m_nBlockSizeBits, (uint64_t) 1) << " to fit window size" << endl;
		numThreads = max(m_nOTs / m_nBlockSizeBits, (uint64_t) 1);
	}

	//The total number of OTs that is performed has to be a multiple of numThreads*Z_REGISTER_BITS
	uint64_t wd_size_bits = m_nBlockSizeBits;//1 << (ceil_log2(m_nBaseOTs));
	uint64_t internal_numOTs = PadToMultiple(ceil_divide(m_nOTs, numThreads), wd_size_bits);

	//Create temporary result buf to which the threads write their temporary masks
	//m_vTempOTMasks.Create(internal_numOTs * numThreads * m_nBitLength);

	//sndthread->Start();
	//rcvthread->Start();

	vector<OTReceiverThread*> rThreads(numThreads);

	for (uint32_t i = 0; i < numThreads; i++) {
		rThreads[i] = new OTReceiverThread(i, internal_numOTs, this);
		rThreads[i]->Start();
	}

	for (uint32_t i = 0; i < numThreads; i++) {
		rThreads[i]->Wait();
	}


	m_nCounter += m_nOTs;

	for (uint32_t i = 0; i < numThreads; i++)
		delete rThreads[i];


	//if (m_eSndOTFlav == Snd_R_OT || m_eSndOTFlav == Snd_GC_OT) {
	//	m_nRet.Copy(m_vTempOTMasks.GetArr(), 0, ceil_divide(m_nOTs * m_nBitLength, 8));
	//}
	//m_vTempOTMasks.delCBitVector();
#ifdef VERIFY_OT
	//Wait for the signal of the corresponding sender thread
	verifyOT(m_nOTs);
#endif

	return true;
}



void OTExtRec::BuildMatrices(CBitVector& T, CBitVector& SndBuf, uint64_t OT_ptr, uint64_t numblocks,
		OT_AES_KEY_CTX* seedkeyptr) {
	uint8_t* ctr_buf = (uint8_t*) calloc (AES_BYTES, sizeof(uint8_t));
	uint64_t* counter = (uint64_t*) ctr_buf;

	uint64_t wd_size_bytes = m_nBlockSizeBytes;//1 << (ceil_log2(m_nBaseOTs) - 3);
	uint64_t rowbytelen = wd_size_bytes * numblocks;
	uint64_t iters = rowbytelen / AES_BYTES;

	uint8_t* Tptr = T.GetArr();
	uint8_t* sndbufptr = SndBuf.GetArr();
	uint8_t* choiceptr;

	//AES_KEY_CTX* seedptr = m_vBaseOTKeys;
	uint64_t global_OT_ptr = OT_ptr + m_nCounter;
	uint32_t keyoffset = m_nBaseOTs;

#ifdef USE_PIPELINED_AES_NI
	//first prg output written to tptrelse
	intrin_sequential_gen_rnd8(ctr_buf, global_OT_ptr, Tptr, iters, (int) m_nBaseOTs, seedkeyptr);

	*counter = global_OT_ptr;
	//second prg output written to snd buffer
	intrin_sequential_gen_rnd8(ctr_buf, global_OT_ptr, sndbufptr, iters, (int) m_nBaseOTs, seedkeyptr+m_nBaseOTs);
#else
	for (uint32_t k = 0; k < m_nBaseOTs; k++) {
		*counter = global_OT_ptr;
		for (uint32_t b = 0; b < iters; b++, (*counter)++) {
			m_cCrypt->encrypt(seedkeyptr + k, Tptr, ctr_buf, AES_BYTES);
			Tptr += AES_BYTES;

			m_cCrypt->encrypt(seedkeyptr + k + keyoffset, sndbufptr, ctr_buf, AES_BYTES);
			sndbufptr += AES_BYTES;
		}
#ifdef DEBUG_OT_SEED_EXPANSION
		cout << "X0[" << k << "]: " << (hex);
		for(uint64_t i = 0; i < AES_BYTES * iters; i++) {
			cout  << setw(2) << setfill('0') << (uint32_t) (Tptr-AES_BYTES*iters)[i];
		}
		cout << (dec) << " (" << (*counter)-iters << ")" <<endl;
		cout << "X1[" << k << "]: " << (hex);
		for(uint64_t i = 0; i < AES_BYTES * iters; i++) {
			cout  << setw(2) << setfill('0') << (uint32_t) (sndbufptr-AES_BYTES*iters)[i];
		}
		cout << (dec) << " (" << (*counter)-iters << ")" <<endl;
#endif
	}
#endif
	//m_vChoices.PrintHex();
	free(ctr_buf);
}

void OTExtRec::MaskBaseOTs(CBitVector& T, CBitVector& SndBuf, uint64_t OTid, uint64_t numblocks) {
	uint64_t rowbytelen = m_nBlockSizeBytes * numblocks;
	uint64_t choicebytelen = bits_in_bytes(min(numblocks * m_nBlockSizeBits, m_nOTs - OTid));
	uint8_t* choiceptr;// = m_nChoices.GetArr() + ceil_divide(OTid, 8);
	CBitVector tmp;


#ifdef GENERATE_T_EXPLICITELY
	//Some nasty moving to compress the code, this part is only required for benchmarking
	uint32_t blockbytesize = rowbytelen * m_nBaseOTs;
	if(m_eRecOTFlav == Rec_R_OT) {
		tmp.CreateBytes(rowbytelen);
		tmp.Reset();
		tmp.XORBytesReverse(SndBuf.GetArr(), 0, rowbytelen);
		tmp.XORBytesReverse(T.GetArr(), 0, rowbytelen);
		m_vChoices->Copy(tmp.GetArr(), ceil_divide(OTid, 8), choicebytelen);

		SndBuf.SetBytes(SndBuf.GetArr()+rowbytelen, blockbytesize-rowbytelen, blockbytesize-rowbytelen);
		SndBuf.SetBytes(T.GetArr()+rowbytelen, 0, blockbytesize-rowbytelen);
		T.FillRand(blockbytesize << 3, m_cCrypt);
		T.SetBytesToZero(0, rowbytelen);
		SndBuf.XORBytes(T.GetArr()+rowbytelen, 0, blockbytesize-rowbytelen);
		SndBuf.XORBytes(T.GetArr()+rowbytelen, blockbytesize-rowbytelen, blockbytesize-rowbytelen);

		for (uint32_t k = 0; k < m_nBaseOTs-1; k++) {
			SndBuf.XORBytesReverse(m_vChoices->GetArr() + ceil_divide(OTid, 8), blockbytesize +  k * rowbytelen, choicebytelen);
		}
	} else {
		uint32_t blockbytesize = rowbytelen * m_nBaseOTs;
		SndBuf.SetBytes(SndBuf.GetArr(), blockbytesize, blockbytesize);
		SndBuf.SetBytes(T.GetArr(), 0, blockbytesize);
		T.FillRand(blockbytesize << 3, m_cCrypt);
		SndBuf.XORBytes(T.GetArr(), 0, blockbytesize);
		SndBuf.XORBytes(T.GetArr(), blockbytesize, blockbytesize);

		for (uint32_t k = 0; k < m_nBaseOTs; k++) {
			SndBuf.XORBytesReverse(m_vChoices->GetArr() + ceil_divide(OTid, 8), blockbytesize +  k * rowbytelen, choicebytelen);
		}
	}

#else
	tmp.CreateBytes(rowbytelen);
	tmp.Reset();

	if(m_eRecOTFlav == Rec_R_OT) {
		tmp.XORBytesReverse(SndBuf.GetArr(), 0, rowbytelen);
		tmp.XORBytesReverse(T.GetArr(), 0, rowbytelen);

		m_vChoices->Copy(tmp.GetArr(), ceil_divide(OTid, 8), choicebytelen);
	} else {
		tmp.Copy(m_vChoices->GetArr() + ceil_divide(OTid, 8), 0, choicebytelen);
	}
	choiceptr = tmp.GetArr();
	for (uint32_t k = 0; k < m_nBaseOTs; k++) {
		SndBuf.XORBytesReverse(choiceptr, k * rowbytelen, rowbytelen);
	}

	SndBuf.XORBytes(T.GetArr(), 0, rowbytelen * m_nBaseOTs);
	tmp.delCBitVector();
#endif
	//cout << "SB: ";
	//SndBuf.PrintHex(0, 32);
}


void OTExtRec::SendMasks(CBitVector Sndbuf, channel* chan, uint64_t OTid, uint64_t processedOTs, uint64_t rec_r_ot_startpos) {
	uint8_t* bufptr = Sndbuf.GetArr();
#ifdef GENERATE_T_EXPLICITELY
	uint64_t nSize = 2 * bits_in_bytes(m_nBaseOTs * processedOTs);
	if(m_eRecOTFlav == Rec_R_OT) {
		nSize = 2 * bits_in_bytes((m_nBaseOTs-1) * processedOTs);
	}
#else
	uint64_t nSize = bits_in_bytes(m_nBaseOTs * processedOTs);
	if(m_eRecOTFlav == Rec_R_OT) {
		nSize = bits_in_bytes((m_nBaseOTs-rec_r_ot_startpos) * processedOTs);
		bufptr = Sndbuf.GetArr() + rec_r_ot_startpos * ceil_divide(processedOTs, 8);
	}
#endif
	chan->send_id_len(bufptr, nSize, OTid, processedOTs);
}


void OTExtRec::HashValues(CBitVector* T, CBitVector* seedbuf, CBitVector* maskbuf, uint64_t OT_ptr, uint64_t OT_len, uint64_t** mat_mul) {
	//uint32_t wd_size_bytes = m_nBlockSizeBytes;//(1 << ((ceil_log2(m_nBaseOTs)) - 3));
	uint32_t rowbytelen = bits_in_bytes(m_nBaseOTs);
	uint32_t hashinbytelen = rowbytelen + sizeof(uint64_t);
	uint32_t aes_key_bytes = m_cCrypt->get_aes_key_bytes();

	uint8_t* Tptr = T->GetArr();
	uint8_t* bufptr = seedbuf->GetArr();

	uint8_t* inbuf = (uint8_t*) calloc(hashinbytelen, 1);
	uint8_t* resbuf = (uint8_t*) calloc(m_cCrypt->get_hash_bytes(), 1);
	uint8_t* hash_buf = (uint8_t*) calloc(m_cCrypt->get_hash_bytes(), 1);


	uint64_t global_OT_ptr = OT_ptr + m_nCounter;
	if(m_eSndOTFlav != Snd_GC_OT) {
		for (uint64_t i = 0; i < OT_len; i++, Tptr += m_nBlockSizeBytes, bufptr += aes_key_bytes, global_OT_ptr++) {
#ifdef DEBUG_OT_HASH_IN
			cout << "Hash-In for i = " << global_OT_ptr << ": " << (hex);
			for(uint32_t p = 0; p < rowbytelen; p++)
				cout << setw(2) << setfill('0') << (uint32_t) Tptr[p];
			cout << (dec) << " (" << (uint32_t) m_vChoices->GetBitNoMask(OT_ptr+i) << ")" << endl;

#endif

#ifdef FIXED_KEY_AES_HASHING
			FixedKeyHashing(m_kCRFKey, bufptr, Tptr, hash_buf, i, ceil_divide(m_nBaseOTs, 8), m_cCrypt);
#else
			memcpy(inbuf, &global_OT_ptr, sizeof(uint64_t));
			memcpy(inbuf+sizeof(uint64_t), Tptr, rowbytelen);
			m_cCrypt->hash_buf(resbuf, aes_key_bytes, inbuf, hashinbytelen, hash_buf);
			memcpy(bufptr, resbuf, aes_key_bytes);
#endif


#ifdef DEBUG_OT_HASH_OUT
			cout << "Hash-Out for i = " << global_OT_ptr << ": " << (hex);
			for(uint32_t p = 0; p < aes_key_bytes; p++)
				cout << setw(2) << setfill('0') << (uint32_t) bufptr[p];
			cout << (dec) << " (" << (uint32_t) m_vChoices->GetBitNoMask(OT_ptr+i) << ")" << endl;
#endif
		}
#ifndef HIGH_SPEED_ROT_LT
		m_fMaskFct->expandMask(maskbuf, seedbuf->GetArr(), 0, OT_len, m_nBitLength, m_cCrypt);

#endif
	} else {
		uint64_t* tmpbuf = (uint64_t*) calloc(PadToMultiple(bits_in_bytes(m_nBitLength), sizeof(uint64_t)), 1);
		uint8_t* tmpbufb = (uint8_t*) calloc(bits_in_bytes(m_nBitLength), 1);

		for(uint64_t i = 0; i < OT_len; i++, Tptr += m_nBlockSizeBytes) {
			BitMatrixMultiplication(tmpbufb, bits_in_bytes(m_nBitLength), Tptr, m_nBaseOTs, mat_mul, tmpbuf);
			//m_vTempOTMasks.SetBits(tmpbufb, (uint64_t) (OT_ptr + i) * m_nBitLength, m_nBitLength);
			maskbuf->SetBits(tmpbufb, i * m_nBitLength, m_nBitLength);
		}
		free(tmpbuf);
		free(tmpbufb);
	}


	free(resbuf);
	free(inbuf);
	free(hash_buf);
}

void OTExtRec::SetOutput(CBitVector* maskbuf, uint64_t otid, uint64_t otlen, queue<mask_block*>* mask_queue,
		channel* chan) {
	uint32_t remots = min(otlen, m_nOTs - otid);

	if (m_eSndOTFlav == Snd_R_OT || m_eSndOTFlav == Snd_GC_OT) {
		CBitVector dummy;//is not used for random OT or GC_OT
		m_fMaskFct->UnMask(otid, remots, m_vChoices, m_vRet, &dummy, maskbuf, m_eSndOTFlav);
	} else {
		mask_block* tmpblock = (mask_block*) malloc(sizeof(mask_block));
		tmpblock->startotid = otid;
		tmpblock->otlen = remots;
		tmpblock->buf = new CBitVector();
		tmpblock->buf->Copy(maskbuf->GetArr(), 0, maskbuf->GetSize());
		//cout << "Creating new tmpblock with startotid = " << otid << " and otlen = " << remots << endl;

		mask_queue->push(tmpblock);
		if(chan->data_available()) {
			ReceiveAndUnMask(chan, mask_queue);
		}
	}
}


void OTExtRec::ReceiveAndUnMask(channel* chan, queue<mask_block*>* mask_queue) {


	uint64_t startotid, otlen, buflen;
	uint8_t *tmpbuf, *buf;
	CBitVector vRcv;
	mask_block* tmpblock;

	while(chan->data_available() && !(mask_queue->empty())) {
		tmpblock = mask_queue->front();
		//Get values and unmask
		buf = chan->blocking_receive_id_len(&tmpbuf, &startotid, &otlen);//chan->blocking_receive();//rcvqueue->front();

		if(startotid != tmpblock->startotid || otlen != tmpblock->otlen) {
			cout << "Startotid = " << startotid << " vs. " << tmpblock->startotid << endl;
			cout << "OTlen = " << otlen << " vs. " << tmpblock->otlen << endl;
		}
		assert(startotid == tmpblock->startotid);
		//cout << " oten = " << otlen << ", tmpblock otlen = " << tmpblock.otlen << endl;
		assert(otlen == tmpblock->otlen);

		buflen = ceil_divide(otlen * m_nBitLength, 8);
		if (m_eSndOTFlav == Snd_OT)
			buflen = buflen * m_nSndVals;
		vRcv.AttachBuf(tmpbuf, buflen);

		uint32_t remots = min(otlen, m_nOTs - startotid);
		m_fMaskFct->UnMask(startotid, remots, m_vChoices, m_vRet, &vRcv, tmpblock->buf, m_eSndOTFlav);
		mask_queue->pop();
		tmpblock->buf->delCBitVector();
		free(buf);
		vRcv.AttachBuf(buf, 0);
		//cout <<  "Start: " << startotid << " data available? " << (uint32_t) chan->data_available() <<
		//		", queue_empty? "<< (uint32_t) mask_queue->empty() << endl;
	}
}

void OTExtRec::ReceiveAndXORCorRobVector(CBitVector& T, uint64_t OT_len, channel* chan) {
	if(m_bUseMinEntCorRob) {
		uint8_t* rndvec = chan->blocking_receive();
		uint64_t len = bits_in_bytes(m_nBaseOTs * OT_len);
		T.XORBytes(rndvec, len);
		free(rndvec);
	}
}


BOOL OTExtRec::verifyOT(uint64_t NumOTs) {
	cout << "Verifying 1oo" << m_nSndVals << " OT" << endl;
	uint32_t nsndvals = m_nSndVals;
	uint32_t choicecodebits = ceil_log2(m_nSndVals);

	CBitVector* vRcvX = new CBitVector[nsndvals];//(CBitVector*) malloc(sizeof(CBitVector)*m_nSndVals);
	for(uint64_t i = 0; i < nsndvals; i++) {
		vRcvX[i].Create(0);
	}
	//vRcvX[1].Create(0);
	CBitVector *Xc, *Xn;
	uint64_t processedOTBlocks, otlen, otstart;
	uint32_t bytelen = ceil_divide(m_nBitLength, 8);
	uint8_t* tempXc = (uint8_t*) malloc(bytelen);
	//uint8_t* tempXn = (uint8_t*) malloc(bytelen);
	//uint8_t** tmpXn = (uint8_t**) malloc(nsndvals-1);
	uint8_t* tempRet = (uint8_t*) malloc(bytelen);
	uint8_t** buf = (uint8_t**) malloc(sizeof(uint8_t*) * nsndvals);
	channel* chan = new channel(OT_ADMIN_CHANNEL, m_cRcvThread, m_cSndThread);
	uint8_t *tmpbuf;
	BYTE resp;
	uint64_t tmpchoice;
	uint64_t tmpret64;

	/*for(uint32_t i = 0; i < nsndvals-1; i++) {
		tmpXn[i] = (uint8_t*) malloc(bytelen);
	}*/

	for (uint64_t i = 0; i < NumOTs;) {
		for(uint64_t j = 0; j < nsndvals; j++) {
			buf[j] = chan->blocking_receive_id_len(&tmpbuf, &otstart, &otlen);
			vRcvX[j].AttachBuf(tmpbuf, bits_in_bytes(otlen * m_nBitLength));
		}

		for (uint64_t j = 0; j < otlen && i < NumOTs; j++, i++) {
			tmpchoice = m_vChoices->Get<uint32_t>(i * choicecodebits, choicecodebits);

			Xc = &vRcvX[tmpchoice];
			Xc->GetBits(tempXc, j * m_nBitLength, m_nBitLength);

			//Xn->GetBits(tempXn, j * m_nBitLength, m_nBitLength);
			m_vRet->GetBits(tempRet, i * m_nBitLength, m_nBitLength);
			for (uint64_t k = 0; k < bytelen; k++) {
				if (tempXc[k] != tempRet[k]) {
					cout << "Error at position i = " << i << ", k = " << k << ", with X" << (hex) << tmpchoice <<
							" = " << (uint32_t) tempXc[k] << " and res = " << (uint32_t) tempRet[k] << ". All values: " <<endl;
					for(uint64_t k = 0, ctr=0; k < nsndvals; k++) {
						tmpret64 = vRcvX[k].Get<uint64_t>(i * m_nBitLength, m_nBitLength);
						cout << "X" << k << ": " << (hex) << tmpret64 << (dec) << endl;
					}

					resp = 0x00;
					chan->send(&resp, 1);

					chan->synchronize_end();
					return false;
				}
			}
		}

		resp = 0x01;
		chan->send(&resp, (uint64_t) 1);

		for(uint64_t j = 0; j < nsndvals; j++) {
			free(buf[j]);
		}
	}

	cout << "OT Verification successful" << endl;

	chan->synchronize_end();
	//cout << "synchronized done" << endl;

	delete chan;
	free(tempXc);
	//free(tempXn);
	free(tempRet);
	free(buf);

	//delete vRcvX;
	return true;
}



void OTExtRec::ComputePKBaseOTs() {
	channel* chan = new channel(OT_ADMIN_CHANNEL, m_cRcvThread, m_cSndThread);
	uint32_t nsndvals = 2;
	uint8_t* pBuf = (uint8_t*) malloc(m_cCrypt->get_hash_bytes() * m_nBaseOTs * nsndvals);
	uint8_t* keyBuf = (uint8_t*) malloc(m_cCrypt->get_aes_key_bytes() * m_nBaseOTs * nsndvals);

	OT_AES_KEY_CTX* tmpkeybuf = (OT_AES_KEY_CTX*) malloc(sizeof(OT_AES_KEY_CTX) * m_nBaseOTs * nsndvals);

	timeval np_begin, np_end;
	gettimeofday(&np_begin, NULL);
	m_cBaseOT->Sender(nsndvals, m_nBaseOTs, chan, pBuf);
	gettimeofday(&np_end, NULL);

#ifdef OTTiming
	printf("Time for performing the base-OTs: %f seconds\n", getMillies(np_begin, np_end));
#endif

	//Key expansion
	uint8_t* pBufIdx = pBuf;
	for(uint32_t i=0; i<m_nBaseOTs; i++ )
	{
		memcpy(keyBuf + i * m_cCrypt->get_aes_key_bytes(), pBufIdx, m_cCrypt->get_aes_key_bytes());
		pBufIdx += m_cCrypt->get_hash_bytes();
		memcpy(keyBuf + i * m_cCrypt->get_aes_key_bytes() + m_nBaseOTs * m_cCrypt->get_aes_key_bytes(), pBufIdx, m_cCrypt->get_aes_key_bytes());
		pBufIdx += m_cCrypt->get_hash_bytes();
	}

	free(pBuf);

	InitPRFKeys(tmpkeybuf, keyBuf, m_nBaseOTs * nsndvals);
	m_tBaseOTKeys.push_back(tmpkeybuf);

	free(keyBuf);
	chan->synchronize_end();

	delete(chan);
}
