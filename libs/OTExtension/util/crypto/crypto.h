/**
 \file 		crypto.h
 \author 	michael.zohner@ec-spride.de
 \copyright __________________
 \brief		Crypto primitive class
 */

#ifndef CRYPTO_H_
#define CRYPTO_H_

#include <openssl/evp.h>
#include <openssl/sha.h>
#include <fstream>
#include <sys/stat.h>
#include <fcntl.h>

#include "../typedefs.h"
#include "../constants.h"
#include "gmp-pk-crypto.h"
#include "ecc-pk-crypto.h"
#include "../socket.h"

#include "TedKrovetzAesNiWrapperC.h"
#include "intrin_sequential_enc8.h"

const uint8_t ZERO_IV[AES_BYTES] = { 0 };

const uint8_t const_seed[] = { 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF };

enum bc_mode {
	ECB, CBC
};

typedef EVP_CIPHER_CTX AES_KEY_CTX;

/* Predefined security levels,
 * ST (SHORTTERM) = 1024/160/163 bit public key, 80 bit private key
 * MT (MEDIUMTERM) = 2048/192/233 bit public key, 112 bit private key
 * LT (LONGTERM) = 3072/256/283 bit public key, 128 bit private key
 * XLT (EXTRA LONGTERM) = 7680/384/409 bit public key, 192 bit private key
 * XXLT (EXTRA EXTRA LONGTERM) = 15360/512/571 bit public key, 256 bit private key
 */

struct prf_state_ctx {
	AES_KEY_CTX aes_key;
	uint64_t* ctr;
};

//TODO: not thread-safe when multiple threads generate random data using the same seed
class crypto {

public:

	crypto(uint32_t symsecbits, uint8_t* seed);
	crypto(uint32_t symsecbits);
	~crypto();

	//Randomness generation routines
	void gen_rnd(uint8_t* resbuf, uint32_t numbytes);
	void gen_rnd_from_seed(uint8_t* resbuf, uint32_t resbytes, uint8_t* seed);
	//void gen_rnd(prf_state_ctx* prf_state, uint8_t* resbuf, uint32_t nbytes);
	void gen_rnd_uniform(uint32_t* res, uint32_t mod);
	void gen_rnd_perm(uint32_t* perm, uint32_t neles);

	//Encryption routines
	void encrypt(uint8_t* resbuf, uint8_t* inbuf, uint32_t ninbytes);
	void decrypt(uint8_t* resbuf, uint8_t* inbuf, uint32_t ninbytes);

	//Hash routines
	void hash(uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes);
	void hash_buf(uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes, uint8_t* buf);
	void hash_non_threadsafe(uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes);
	void hash_ctr(uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes, uint64_t ctr);
	void fixed_key_aes_hash(AES_KEY_CTX* aes_key, uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes);
	void fixed_key_aes_hash_ctr(uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes);

	//Key seed routines
	void seed_aes_hash(uint8_t* seed, bc_mode mode = ECB, const uint8_t* iv = ZERO_IV);
	void seed_aes_enc(uint8_t* seed, bc_mode mode = ECB, const uint8_t* iv = ZERO_IV);

	//External encryption routines
	void init_aes_key(AES_KEY_CTX* aes_key, uint8_t* seed, bc_mode mode = ECB, const uint8_t* iv = ZERO_IV);
	void init_aes_key(AES_KEY_CTX* aes_key, uint32_t symbits, uint8_t* seed, bc_mode mode = ECB, const uint8_t* iv = ZERO_IV);
	void clean_aes_key(AES_KEY_CTX* aeskey);
	uint32_t get_aes_key_bytes();
	void encrypt(AES_KEY_CTX* enc_key, uint8_t* resbuf, uint8_t* inbuf, uint32_t ninbytes);
	void decrypt(AES_KEY_CTX* dec_key, uint8_t* resbuf, uint8_t* inbuf, uint32_t ninbytes);

	pk_crypto* gen_field(field_type ftype);

	seclvl get_seclvl() {
		return secparam;
	}
	;
	uint32_t get_hash_bytes();

	void gen_common_seed(prf_state_ctx* aes_key, CSocket& sock);

	void init_prf_state(prf_state_ctx* prf_state, uint8_t* seed);
	void free_prf_state(prf_state_ctx* prf_state);
private:
	void seed_aes_key(AES_KEY_CTX* aeskey, uint8_t* seed, bc_mode mode = ECB, const uint8_t* iv = ZERO_IV, bool encrypt = true);
	void seed_aes_key(AES_KEY_CTX* aeskey, uint32_t symseclvl, uint8_t* seed, bc_mode mode = ECB, const uint8_t* iv = ZERO_IV, bool encrypt = true);
	void init(uint32_t symsecbits, uint8_t* seed);


	AES_KEY_CTX aes_hash_key;
	AES_KEY_CTX aes_enc_key;
	AES_KEY_CTX aes_dec_key;
	prf_state_ctx global_prf_state;

	seclvl secparam;
	uint8_t* aes_hash_in_buf;
	uint8_t* aes_hash_out_buf;
	uint8_t* aes_hash_buf_y1;
	uint8_t* aes_hash_buf_y2;

	uint8_t* sha_hash_buf;

	void (*hash_routine)(uint8_t*, uint32_t, uint8_t*, uint32_t, uint8_t*);
};

//Some functions that should be useable without the class
void sha1_hash(uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes, uint8_t* hash_buf);
void sha256_hash(uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes, uint8_t* hash_buf);
void sha512_hash(uint8_t* resbuf, uint32_t noutbytes, uint8_t* inbuf, uint32_t ninbytes, uint8_t* hash_buf);
void gen_secure_random(uint8_t* dest, uint32_t nbytes);
void gen_rnd_bytes(prf_state_ctx* prf_state, uint8_t* resbuf, uint32_t nbytes);

seclvl get_sec_lvl(uint32_t symsecbits); //TODO pick a more elegant name (see crypto->get_seclvl())

//static const uint32_t m_nCodeWordBits = 256;
//static const uint32_t m_nCodeWordBytes = m_nCodeWordBits / 8;

#endif /* CRYPTO_H_ */
