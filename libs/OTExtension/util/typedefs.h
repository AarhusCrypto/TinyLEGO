#ifndef __TYPEDEFS_H__BY_SGCHOI
#define __TYPEDEFS_H__BY_SGCHOI

#include <sys/time.h>
#include <assert.h>
#include <float.h>



#define two_pow(e) (((uint64_t) 1) << (e))

static int ceil_log2(int bits) {
	if (bits == 1)
		return 1;
	int targetlevel = 0, bitstemp = bits;
	while (bitstemp >>= 1)
		++targetlevel;
	return targetlevel + ((1 << targetlevel) < bits);
}

static int floor_log2(int bits) {
	if (bits == 1)
		return 1;
	int targetlevel = 0;
	while (bits >>= 1)
		++targetlevel;
	return targetlevel;
}
#define pad_to_power_of_two(e) ( ((uint64_t) 1) << (ceil_log2(e)) )

static double getMillies(timeval timestart, timeval timeend) {
	long time1 = (timestart.tv_sec * 1000000) + (timestart.tv_usec);
	long time2 = (timeend.tv_sec * 1000000) + (timeend.tv_usec);

	return (double) (time2 - time1) / 1000;
}


enum e_role {
	SERVER, CLIENT, ALL
};

typedef struct SECURITYLEVELS {
	int statbits;
	int symbits;
	int ifcbits;
	int eccpfbits;
	int ecckcbits;
} seclvl;

typedef int BOOL;
typedef long LONG;

typedef unsigned char BYTE;
typedef unsigned short USHORT;
typedef unsigned int UINT;
typedef unsigned long ULONG;
typedef BYTE UINT8_T;
typedef USHORT UINT16_T;
typedef UINT UINT32_T;
typedef unsigned long long UINT64_T;
typedef long long SINT64_T;

typedef ULONG DWORD;
typedef UINT64_T UGATE_T;
typedef UINT64_T REGISTER_SIZE;

#define GATE_T_BITS (sizeof(UGATE_T) * 8)

typedef REGISTER_SIZE REGSIZE;
#define LOG2_REGISTER_SIZE		ceil_log2(sizeof(REGISTER_SIZE) << 3)

#define FILL_BYTES				AES_BYTES
#define FILL_BITS				AES_BITS

#define OT_WINDOW_SIZE		(AES_BITS*4)
#define OT_WINDOW_SIZE_BYTES	(AES_BYTES*4)

#define MAX_REPLY_BITS			65536 //at most 2^16 bits may be sent in one go

#define RETRY_CONNECT		1000
#define CONNECT_TIMEO_MILISEC	10000

#define SNDVALS 2

#define OTEXT_BLOCK_SIZE_BITS	AES_BITS
#define OTEXT_BLOCK_SIZE_BYTES	AES_BYTES

#define VECTOR_INTERNAL_SIZE 8

#define			SERVER_ID	0
#define			CLIENT_ID	1

#define MAX_INT (~0)
#if (MAX_INT == 0xFFFFFFFF)
#define MACHINE_SIZE_32
#elif (MAX_INT == 0xFFFFFFFFFFFFFFFF)
#define MACHINE_SIZE_64
#else
#define MACHINE_SIZE_16
#endif

template<class T>
T rem(T a, T b) {
	return ((a) > 0) ? (a) % (b) : (a) % (b) + ((b) > 0 ? (b) : (b) * -1);
}
template<class T>
T sub(T a, T b, T m) {
	return ((b) > (a)) ? (a) + (m) - (b) : (a) - (b);
}
#ifndef FALSE
#define FALSE			0
#endif
#ifndef TRUE
#define TRUE			1
#endif
#define ZERO_BYTE		0
#define MAX_BYTE		0xFF
#define MAX_UINT		0xFFFFFFFF

#ifdef WIN32
#include <WinSock2.h>
#include <windows.h>

typedef unsigned short USHORT;
typedef int socklen_t;
#pragma comment(lib, "wsock32.lib")

#define SleepMiliSec(x)			Sleep(x)

#else //WIN32

#include <sys/types.h>       
#include <sys/socket.h>      
#include <netdb.h>           
#include <arpa/inet.h>       
#include <unistd.h>          
#include <netinet/in.h>   
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <netinet/tcp.h>
#include <queue>

typedef int SOCKET;
#define INVALID_SOCKET -1

#define SleepMiliSec(x)			usleep((x)<<10)
#endif// WIN32

#define ceil_divide(x, y)			(( ((x) + (y)-1)/(y)))
#define bits_in_bytes(bits) (ceil_divide((bits), 8))


#define PadToRegisterSize(x) 		(PadToMultiple(x, OTEXT_BLOCK_SIZE_BITS))
#define PadToMultiple(x, y) 		( ceil_divide(x, y) * (y))

#include <cstring>
#include <string>  
#include <vector> 
#include <iostream>

using namespace std;

#endif //__TYPEDEFS_H__



