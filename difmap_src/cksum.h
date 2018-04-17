#ifndef cksum_h
#define cksum_h

/*
 * The CheckSum type holds a table of CRC byte checksums to facilitate up
 * the process of calculating the checksum of an array of bytes.
 */
typedef struct CheckSum CheckSum;

CheckSum *new_CheckSum(void);
CheckSum *del_CheckSum(CheckSum *cs);

/*
 * Compute a 32-bit cyclic redundancy checksum of an object of given
 * size, using the ethernet key.
 */
unsigned long checksum_of_object(CheckSum *cs, void *obj, size_t nbyte);

#endif
