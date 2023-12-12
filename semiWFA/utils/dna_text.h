#ifndef DNA_TEXT_H_
#define DNA_TEXT_H_
#include "utils/commons.h"
#define DNA_RANGE           4
#define DNA_EXTENDED_RANGE  5
#define DNA_RANGE_BITS 2
#define DNA_CHAR_A 'A'
#define DNA_CHAR_C 'C'
#define DNA_CHAR_G 'G'
#define DNA_CHAR_T 'T'
#define DNA_CHAR_N 'N'
#define ENC_DNA_CHAR_A 0
#define ENC_DNA_CHAR_C 1
#define ENC_DNA_CHAR_G 2
#define ENC_DNA_CHAR_T 3
#define ENC_DNA_CHAR_N 4
extern const uint8_t dna_encode_table[256];
extern const char dna_decode_table[DNA_EXTENDED_RANGE];
#define dna_encode(character)     (dna_encode_table[(int)(character)])
#define dna_decode(enc_char)      (dna_decode_table[(int)(enc_char)])
#endif 
