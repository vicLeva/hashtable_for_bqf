#include <cstdint>
#include <string>
#include <iostream>
#include "additionnal_methods.hpp"

void print_bits(uint64_t x) {
  std::bitset<MEM_UNIT> bits(x);
  std::cout << bits << std::endl;
}

uint64_t mask_right(uint64_t numbits){
    uint64_t mask = -(numbits >= MEM_UNIT) | ((1ULL << numbits) - 1ULL);
    return mask;
}

uint64_t encode(std::string kmer){
    uint64_t encoded = 0;
    for(char& c : kmer) {
        if (c=='T'){
            encoded <<= 2;
            encoded |= 3;
        }
        else if (c=='G'){
            encoded <<= 2;
            encoded |= 2;
        }
        else if (c=='C'){
            encoded <<= 2;
            encoded |= 1;
        }
        else{ //T is 00 so with reverse complementarity we won't get 0000000000000 as input for xorshift
            encoded <<= 2;
        }
    }

    return encoded;
}

std::string decode(uint64_t revhash, uint64_t k){
    std::string kmer;

    for (size_t i=0; i<k; i++){
        switch(revhash & mask_right(2)){
            case 3:
                kmer = 'T' + kmer;
                break;
            case 2:
                kmer = 'G' + kmer;
                break;
            case 1:
                kmer = 'C' + kmer;
                break;
            default:
                kmer = 'A' + kmer;
        }
        revhash >>= 2;
    }

    return kmer;
}


uint64_t nucl_encode(char nucl){
  //Returns the binary encoding of a nucleotide
  //different from encode() function because this is for query, so we have to check for canonical smer
  //and this encoding allows fast comparison (<) of lexico order 
  switch (nucl){
    case 'A':
      return 0;
    case 'T':
      return 3;
    case 'C':
      return 1;
    case 'G':
      return 2;
    default :
        std::cout << "non nucl : " << nucl << "\n";
        throw std::invalid_argument( "received non nucleotidic value" );
  }
}

uint64_t flip(uint64_t encoding, size_t bitsize){
  //Used so Ts are 0s so we dont get 0000.. in xorshift (through revcompl)
    return ~encoding << (64-bitsize) >> (64-bitsize);
}

uint64_t revcomp64 (const uint64_t v, size_t bitsize){
  return (((uint64_t)rev_table[v & 0xff] << 56) | 
    ((uint64_t)rev_table[(v >> 8) & 0xff] << 48) | 
    ((uint64_t)rev_table[(v >> 16) & 0xff] << 40) | 
    ((uint64_t)rev_table[(v >> 24) & 0xff] << 32) | 
    ((uint64_t)rev_table[(v >> 32) & 0xff] << 24) | 
    ((uint64_t)rev_table[(v >> 40) & 0xff] << 16) |
    ((uint64_t)rev_table[(v >> 48) & 0xff] << 8) |
    ((uint64_t)rev_table[(v >> 56) & 0xff])) >> (64-bitsize);
}

uint64_t canonical(uint64_t smer, size_t size){
    uint64_t revcomp = revcomp64(smer, size);
    if (revcomp < smer) { return revcomp; }
    else { return smer; }
}

std::string canonical(const std::string& smer, size_t s){
  return decode(canonical(encode(smer), 2*s), s);
}