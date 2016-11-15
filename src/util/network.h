#ifndef TINY_UTIL_NETWORK_H_
#define TINY_UTIL_NETWORK_H_

#include "OTExtension/util/socket.h"
#include "OTExtension/util/rcvthread.h"
#include "OTExtension/util/sndthread.h"
#include "OTExtension/util/channel.h"

#define RETRY_CONNECT   1000
#define CONNECT_TIMEO_MILISEC 10000

class Network {
public:
  Network(const char m_nAddr[], USHORT m_nPort);

  void ConnectAndStart();
  void ListenAndStart();
  BOOL Connect();
  BOOL Listen();

  const char* m_nAddr;
  USHORT m_nPort;
  CSocket* m_vSocket;
  SndThread* sndthread;
  RcvThread* rcvthread;
};

#endif /* TINY_UTIL_NETWORK_H_ */