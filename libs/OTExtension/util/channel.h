/*
 * channel.h
 *
 *  Created on: Mar 9, 2015
 *      Author: mzohner
 */

#ifndef CHANNEL_H_
#define CHANNEL_H_

#include "typedefs.h"
#include "socket.h"
#include "rcvthread.h"
#include "sndthread.h"

class channel {
public:
	channel(uint16_t channelid, RcvThread* rcver, SndThread* snder) {
		m_cRcver = rcver;
		m_cSnder = snder;
		m_bChannelID = channelid;
		m_eRcved = new CEvent;
		m_eFin = new CEvent;
		m_qRcvedBlocks = rcver->add_listener(channelid, m_eRcved, m_eFin);
		m_bSndAlive = true;
		m_bRcvAlive = true;
	}

	~channel() {
		if(m_bRcvAlive) {
			m_cRcver->remove_listener(m_bChannelID);
		}

		delete m_eRcved;
		delete m_eFin;
	}

	void send(uint8_t* buf, uint64_t nbytes) {
		assert(m_bSndAlive);
		m_cSnder->add_snd_task(m_bChannelID, nbytes, buf);
	}
	void send_id_len(uint8_t* buf, uint64_t nbytes, uint64_t id, uint64_t len) {
		assert(m_bSndAlive);
		m_cSnder->add_snd_task_start_len(m_bChannelID, nbytes, buf, id, len);
	}

	//buf needs to be freed, data contains the payload
	uint8_t* blocking_receive_id_len(uint8_t** data, uint64_t* id, uint64_t* len) {
		uint8_t* buf = blocking_receive();
		*data = buf;
		*id = *((uint64_t*) *data);
		(*data)  += sizeof(uint64_t);
		*len = *((uint64_t*) *data);
		(*data) += sizeof(uint64_t);

		return buf;
	}

	uint8_t* blocking_receive() {
		assert(m_bRcvAlive);
		while(m_qRcvedBlocks->empty())
			m_eRcved->Wait();
		uint8_t* ret_block = m_qRcvedBlocks->front();
		m_qRcvedBlocks->pop();

		return ret_block;
	}

	bool is_alive() {
		return (!(m_qRcvedBlocks->empty() && m_eFin->IsSet()));
	}

	bool data_available() {
		return !m_qRcvedBlocks->empty();
	}

	void signal_end() {
		m_cSnder->signal_end(m_bChannelID);
		m_bSndAlive = false;
	}

	void wait_for_fin() {
		m_eFin->Wait();
		m_bRcvAlive = false;
	}

	//TODO
	void synchronize() {

	}

	void synchronize_end() {
		if(m_bSndAlive)
			signal_end();
		if(m_bRcvAlive)
			m_cRcver->flush_queue(m_bChannelID);
		if(m_bRcvAlive)
			wait_for_fin();

	}

private:
	RcvThread* m_cRcver;
	SndThread* m_cSnder;
	CEvent* m_eRcved;
	CEvent* m_eFin;
	uint16_t m_bChannelID;
	queue<uint8_t*>* m_qRcvedBlocks;
	bool m_bSndAlive;
	bool m_bRcvAlive;
};


#endif /* CHANNEL_H_ */
