#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace Simp21 {

//-------------------------------------------------------------------------------------------------
// A simple buffered character stream class:

static const int buffer_size = 1048576;

class StreamBuffer {
  FILE *in;
  unsigned char buf[buffer_size];
  int pos;
  int size;

  void assureLookahead() {
    if (pos >= size) {
      pos = 0;
      //             size = gzread(in, buf, sizeof(buf));
      size = fread(buf, 1, sizeof(buf), in);
    }
  }

public:
  explicit StreamBuffer(FILE *i) : in(i), pos(0), size(0) { assureLookahead(); }

  int operator*() const { return (pos >= size) ? EOF : buf[pos]; }
  void operator++() {
    pos++;
    assureLookahead();
  }
  int position() const { return pos; }
};

//-------------------------------------------------------------------------------------------------
// End-of-file detection functions for StreamBuffer and char*:

static inline bool isEof(StreamBuffer &in) { return *in == EOF; }
static inline bool isEof(const char *in) { return *in == '\0'; }

//-------------------------------------------------------------------------------------------------
// Generic parse functions parametrized over the input-stream type.

template <class B> static void skipWhitespace(B &in) {
  while ((*in >= 9 && *in <= 13) || *in == 32) {
    ++in;
  }
}

template <class B> static void skipLine(B &in) {
  for (;;) {
    if (isEof(in)) {
      return;
    }

    if (*in == '\n') {
      ++in;
      return;
    }

    ++in;
  }
}

template <class B> static double parseDouble(B &in) { // only in the form X.XXXXXe-XX
  bool neg = false;
  double accu = 0.0;
  double currentExponent = 1;
  int exponent;
  skipWhitespace(in);

  if (*in == EOF) {
    return 0;
  }

  if (*in == '-') {
    neg = true, ++in;
  } else if (*in == '+') {
    ++in;
  }

  if (*in < '1' || *in > '9') {
    printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  }

  accu = (double)(*in - '0');
  ++in;

  if (*in != '.') {
    printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  }

  ++in; // skip dot
  currentExponent = 0.1;

  while (*in >= '0' && *in <= '9') {
    accu = accu + currentExponent * ((double)(*in - '0')), currentExponent /= 10, ++in;
  }

  if (*in != 'e') {
    printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  }

  ++in;                    // skip dot
  exponent = parseInt(in); // read exponent
  accu *= pow(10, exponent);
  return neg ? -accu : accu;
}

template <class B> static int parseInt(B &in) {
  int val = 0;
  bool neg = false;
  skipWhitespace(in);

  if (*in == '-') {
    neg = true, ++in;
  } else if (*in == '+') {
    ++in;
  }

  if (*in < '0' || *in > '9') {
    fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  }

  while (*in >= '0' && *in <= '9')
    val = val * 10 + (*in - '0'), ++in;

  return neg ? -val : val;
}

// String matching: in case of a match the input iterator will be advanced the corresponding
// number of characters.
template <class B> static bool match(B &in, const char *str) {
  int i;

  for (i = 0; str[i] != '\0'; i++)
    if (in[i] != str[i]) {
      return false;
    }

  in += i;
  return true;
}

// String matching: consumes characters eagerly, but does not require random access iterator.
template <class B> static bool eagerMatch(B &in, const char *str) {
  for (; *str != '\0'; ++str, ++in)
    if (*str != *in) {
      return false;
    }

  return true;
}

} // namespace Simp21
