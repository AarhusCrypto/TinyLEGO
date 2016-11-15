#include "util/network.h"

Network::Network(const char* m_nAddr, USHORT m_nPort) : m_nAddr(m_nAddr) {
  m_vSocket = new CSocket();

  sndthread = new SndThread(m_vSocket);
  rcvthread = new RcvThread(m_vSocket);

  this->m_nPort = m_nPort;
}

void Network::ConnectAndStart() {
  Connect();
  sndthread->Start();
  rcvthread->Start();
}

void Network::ListenAndStart() {
  Listen();
  sndthread->Start();
  rcvthread->Start();
}

BOOL Network::Connect() {
  bool bFail = FALSE;
  uint64_t lTO = CONNECT_TIMEO_MILISEC;

  for (int k = 0; k >= 0 ; k--)
  {
    for ( int i = 0; i < RETRY_CONNECT; i++ )
    {
      if ( !m_vSocket->Socket() )
      {
        printf("Socket failure: ");
        goto connect_failure;
      }

      if ( m_vSocket->Connect( m_nAddr, m_nPort, lTO))
      {
        // send pid when connected
        m_vSocket->Send( &k, sizeof(int) );
        if (k == 0)
        {
          //std::cout << "connected" << std::endl;
          return TRUE;
        }
        else
        {
          break;
        }
        SleepMiliSec(10);
        m_vSocket->Close();
      }
      SleepMiliSec(20);
      if (i + 1 == RETRY_CONNECT)
        goto server_not_available;
    }
  }
server_not_available:
  printf("Server not available: ");
connect_failure:
  return FALSE;
}

BOOL Network::Listen() {
#ifdef TINY_PRINT
  std::cout << "Listening: " << m_nAddr << ":" << m_nPort << ", with size: " << std::endl;
#endif
  if ( !m_vSocket->Socket() )
  {
    goto listen_failure;
  }
  if ( !m_vSocket->Bind(m_nPort, m_nAddr) )
    goto listen_failure;
  if ( !m_vSocket->Listen() )
    goto listen_failure;

  for ( int i = 0; i < 1; i++ ) //twice the actual number, due to double sockets for OT
  {
    CSocket sock;
    //std::cout << "New round! " << std::endl;
    if ( !m_vSocket->Accept(sock) )
    {
      std::cerr << "Error in accept" << std::endl;
      goto listen_failure;
    }

    UINT threadID;
    sock.Receive(&threadID, sizeof(int));

    if ( threadID >= 1)
    {
      sock.Close();
      i--;
      continue;
    }

#ifdef TINY_PRINT
    std::cout <<  " (" << role << ") (" << threadID << ") connection accepted" << std::endl;
#endif
    // locate the socket appropriately
    m_vSocket->AttachFrom(sock);
    sock.Detach();
  }

#ifdef TINY_PRINT
  std::cout << "Listening finished"  << std::endl;
#endif
  return TRUE;

listen_failure:
  std::cout << "Listen failed" << std::endl;
  return FALSE;
}