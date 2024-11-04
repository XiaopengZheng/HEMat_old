#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include "2024.h"
#include <helib/helib.h>
#include <vector>
#include <cmath>

using namespace std;
using namespace NTL;

void matrix_nosquare(int m, int l, int n)
{
  long p;
  p = 202001;
  unsigned long m0 = 101 * 5 * 5 * 5;
  unsigned long r = 1;
  int l_0 = l, m_0 = m, n_0 = n;
  unsigned long bits = 120;
  unsigned long c = 2;
  const int slot_b = 100;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m0
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  std::cout << "  Generating key-switching matrices..." << std::endl;
  HELIB_NTIMER_START(KeyGen);

  secret_key.GenSecKey(); //
  int down = 2;
  int right = 11;
  secret_key.GenKeySWmatrix(1, 2, 0, 0);
  secret_key.GenKeySWmatrix(1, 11, 0, 0);
  secret_key.setKeySwitchMap();

  const helib::PubKey &public_key = secret_key;
  vector<int> factors = {1, 2, 4, 5, 10, 20, 25, 50, 100};
  int D0, D1;
  context.printout();
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(m, l, 10, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M2 is equal to: " << std::endl;
  std::vector<std::vector<long>> M2 = randmat(l, n, 10, p, "M2.txt");
  std::cout << std::endl;
  cout << "Plaintext multiplication: M1*M2 mod " << p << " is equal to: " << endl;
  plain_text_multiplication(M1, M2, m, l, n, p);
  if (l <= m && l <= n)
  {
    for (int i = 0; i < factors.size(); i++)
    {
      if (factors[i] >= l)
      {
        l = factors[i];
        break;
      }
    }
    for (int i = factors.size() - 1; i >= 0; i--)
    {
      if (m <= factors[i] && (factors[i] % l == 0))
      {
        D0 = factors[i];
        m = factors[i];
      }
      if (n <= factors[i] && ((factors[i] % l) == 0))
      {
        D1 = factors[i];
        n = factors[i];
      }
    }
    std::vector<long> LeftRot, UpRot;
    LeftRot.push_back(1);
    UpRot.push_back(1);
    int l0 = slot_b / D0;
    int l1 = slot_b / D1;
    for (int i = 1; i < D1; i++)
    {
      LeftRot.push_back(Modpower(right, l1 * (D1 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(right, l1 * (D1 - i), m0), 0, 0);
    }
    for (int i = 1; i < D0; i++)
    {
      UpRot.push_back(Modpower(down, l0 * (D0 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(down, l0 * (D0 - i), m0), 0, 0);
    }
    secret_key.setKeySwitchMap();
    // for (int i = 1; i < D0; i++)
    // {
    //   RightRot.push_back(Modpower(3, b_size * b_size * i - b_size * i, m));
    //   secret_key.GenKeySWmatrix(1, Modpower(3, b_size * b_size * i - b_size * i, m), 0, 0);
    // }
    helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
    // Encode
    Encode(ptxt1, ptxt2, M1, M2, slot_b, l0, l1);

    //  Encrypt
    helib::Ctxt ctxt1(public_key), ctxt2(public_key);

    public_key.Encrypt(ctxt1, ptxt1);
    public_key.Encrypt(ctxt2, ptxt2);
    HELIB_NTIMER_START(Mult_2023b);
    ctxt1 = Replicate(ctxt1, D1 / l, l, LeftRot);
    ctxt2 = Replicate(ctxt2, D0 / l, l, UpRot);
    helib::Ctxt ctxt_result(public_key);
    ctxt_result = FHEMatMultMain(context, ctxt1, ctxt2, m, l, n, D0, D1, slot_b, l0, l1, LeftRot, UpRot);
    HELIB_NTIMER_STOP(Mult_2023b);
    helib::printNamedTimer(std::cout, "Mult_2023b");
    //  Decrypt:
    helib::Ptxt<helib::BGV> plaintext_result(context);
    secret_key.Decrypt(plaintext_result, ctxt_result);
    printf_ptxt(plaintext_result, D0, D1, slot_b, m_0, n_0, p);
    std::cout << std::endl;
  }
  else if (m >= l && l >= n)
  {
    D0 = m;
    D1 = l;
    std::vector<long> LeftRot, UpRot;
    LeftRot.push_back(1);
    UpRot.push_back(1);
    int l0 = slot_b / D0;
    int l1 = slot_b / D1;
    for (int i = 1; i < D1; i++)
    {
      LeftRot.push_back(Modpower(right, l1 * (D1 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(right, l1 * (D1 - i), m0), 0, 0);
    }
    for (int i = 1; i < D0; i++)
    {
      UpRot.push_back(Modpower(down, l0 * (D0 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(down, l0 * (D0 - i), m0), 0, 0);
    }
    secret_key.setKeySwitchMap();
    helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
    // Encode
    Encode(ptxt1, ptxt2, M1, M2, slot_b, l0, l1);
    // Encrypt

    helib::Ctxt ctxt1(public_key), ctxt2(public_key);
    public_key.Encrypt(ctxt1, ptxt1);
    public_key.Encrypt(ctxt2, ptxt2);
    //
    HELIB_NTIMER_START(Mult_2023b);
    ctxt2 = Replicate(ctxt2, D1 / n, n, LeftRot);
    ctxt2 = Replicate(ctxt2, D0 / l, l, UpRot);
    helib::Ctxt ctxt_result(public_key);
    ctxt_result = FHEMatMultMain(context, ctxt1, ctxt2, m, l, n, D0, D1, slot_b, l0, l1, LeftRot, UpRot);
    ctxt_result = Replicate(ctxt_result, D1 / n, n, LeftRot);
    HELIB_NTIMER_STOP(Mult_2023b);
    helib::printNamedTimer(std::cout, "Mult_2023b");
    // Decrypt
    helib::Ptxt<helib::BGV> plaintext_result(context);
    secret_key.Decrypt(plaintext_result, ctxt_result);
    printf_ptxt(plaintext_result, D0, D1, slot_b, m_0, n_0, p);
  }
  else if (n >= l && l >= m)
  {
    D0 = l;
    D1 = n;
    std::vector<long> LeftRot, UpRot;
    LeftRot.push_back(1);
    UpRot.push_back(1);
    int l0 = slot_b / D0;
    int l1 = slot_b / D1;
    for (int i = 1; i < D1; i++)
    {
      LeftRot.push_back(Modpower(right, l1 * (D1 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(right, l1 * (D1 - i), m0), 0, 0);
    }
    for (int i = 1; i < D0; i++)
    {
      UpRot.push_back(Modpower(down, l0 * (D0 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(down, l0 * (D0 - i), m0), 0, 0);
    }
    secret_key.setKeySwitchMap();
    helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
    // Encode
    Encode(ptxt1, ptxt2, M1, M2, slot_b, l0, l1);
    // Encrypt

    helib::Ctxt ctxt1(public_key), ctxt2(public_key);
    public_key.Encrypt(ctxt1, ptxt1);
    public_key.Encrypt(ctxt2, ptxt2);
    //
    HELIB_NTIMER_START(Mult_2023b);
    ctxt1 = Replicate(ctxt1, D1 / l, l, LeftRot);
    ctxt1 = Replicate(ctxt1, D0 / m, m, UpRot);
    helib::Ctxt ctxt_result(public_key);
    ctxt_result = FHEMatMultMain(context, ctxt1, ctxt2, m, l, n, D0, D1, slot_b, l0, l1, LeftRot, UpRot);
    ctxt_result = Replicate(ctxt_result, D0 / m, m, UpRot);
    HELIB_NTIMER_STOP(Mult_2023b);
    helib::printNamedTimer(std::cout, "Mult_2023b");
    // Decrypt
    helib::Ptxt<helib::BGV> plaintext_result(context);
    secret_key.Decrypt(plaintext_result, ctxt_result);
    printf_ptxt(plaintext_result, D0, D1, slot_b, m_0, n_0, p);
  }
  else if (l >= m && l >= n)
  {
    D0 = l;
    D1 = l;
    std::vector<long> LeftRot, UpRot;
    LeftRot.push_back(1);
    UpRot.push_back(1);
    int l0 = slot_b / D0;
    int l1 = slot_b / D1;
    for (int i = 1; i < D1; i++)
    {
      LeftRot.push_back(Modpower(right, l1 * (D1 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(right, l1 * (D1 - i), m0), 0, 0);
    }
    for (int i = 1; i < D0; i++)
    {
      UpRot.push_back(Modpower(down, l0 * (D0 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(down, l0 * (D0 - i), m0), 0, 0);
    }
    secret_key.setKeySwitchMap();
    helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
    // Encode
    Encode(ptxt1, ptxt2, M1, M2, slot_b, l0, l1);
    // Encrypt

    helib::Ctxt ctxt1(public_key), ctxt2(public_key);
    public_key.Encrypt(ctxt1, ptxt1);
    public_key.Encrypt(ctxt2, ptxt2);
    //
    HELIB_NTIMER_START(Mult_2023b);
    ctxt1 = Replicate(ctxt1, D0 / m, m, UpRot);
    ctxt2 = Replicate(ctxt2, D1 / n, n, LeftRot);
    helib::Ctxt ctxt_result(public_key);
    ctxt_result = FHEMatMultMain(context, ctxt1, ctxt2, m, l, n, D0, D1, slot_b, l0, l1, LeftRot, UpRot);
    // ctxt_result = Replicate(ctxt_result, D1 / n, n, LeftRot);
    ctxt_result = Replicate(ctxt_result, D0 / m, m, UpRot);
    HELIB_NTIMER_STOP(Mult_2023b);
    helib::printNamedTimer(std::cout, "Mult_2023b");
    // Decrypt
    helib::Ptxt<helib::BGV> plaintext_result(context);
    secret_key.Decrypt(plaintext_result, ctxt_result);
    printf_ptxt(plaintext_result, D0, D1, slot_b, m_0, n_0, p);
  }
}

void example(int m, int l, int n)
{
  long p = 3221319251;
  unsigned long m0 = 101 * 5 * 5 * 5 * 3;
  unsigned long r = 1;
  int l_0 = l, m_0 = m, n_0 = n;
  unsigned long bits = 530;
  unsigned long c = 2;
  const int slot_b = 100;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m0
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  std::cout << "  Generating key-switching matrices..." << std::endl;
  HELIB_NTIMER_START(KeyGen);
  // HeapProfilerStart("test"); //开始监测

  secret_key.GenSecKey(); //
  int down = 2;
  int right = 11;
  secret_key.GenKeySWmatrix(1, 2, 0, 0);
  secret_key.GenKeySWmatrix(1, 11, 0, 0);
  secret_key.setKeySwitchMap();

  const helib::PubKey &public_key = secret_key;
  vector<int> factors = {1, 2, 4, 5, 10, 20, 25, 50, 100};
  int D0, D1;
  context.printout();
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(m, l, 10, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M2 is equal to: " << std::endl;
  std::vector<std::vector<long>> M2 = randmat(l, n, 10, p, "M2.txt");
  std::cout << std::endl;
  cout << "Plaintext multiplication: M1*M2*M2*M2 mod " << p << " is equal to: " << endl;
  std::vector<std::vector<long>> M;
  M = plain_text_multiplication_output(M1, M2, m, l, n, p);
  M = plain_text_multiplication_output(M, M2, m, l, n, p);
  M = plain_text_multiplication_output(M, M2, m, l, n, p);
  M = plain_text_multiplication_output(M, M2, m, l, n, p);
  printf_matrix(M, m, n, log10(p) + 2, p);
  D0 = l;
  D1 = n;
  std::vector<long> LeftRot, UpRot;
  LeftRot.push_back(1);
  UpRot.push_back(1);
  int l0 = slot_b / D0;
  int l1 = slot_b / D1;
  for (int i = 1; i < D1; i++)
  {
    LeftRot.push_back(Modpower(right, l1 * (D1 - i), m0));
    secret_key.GenKeySWmatrix(1, Modpower(right, l1 * (D1 - i), m0), 0, 0);
  }
  for (int i = 1; i < D0; i++)
  {
    UpRot.push_back(Modpower(down, l0 * (D0 - i), m0));
    secret_key.GenKeySWmatrix(1, Modpower(down, l0 * (D0 - i), m0), 0, 0);
  }
  secret_key.setKeySwitchMap();
  helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
  // Encode
  Encode(ptxt1, ptxt2, M1, M2, slot_b, l0, l1);
  // Encrypt

  helib::Ctxt ctxt1(public_key), ctxt2(public_key);
  public_key.Encrypt(ctxt1, ptxt1);
  public_key.Encrypt(ctxt2, ptxt2);
  //
  HELIB_NTIMER_START(Mult_2023b);
  ctxt1 = Replicate(ctxt1, D1 / l, l, LeftRot);
  ctxt1 = Replicate(ctxt1, D0 / m, m, UpRot);
  helib::Ctxt ctxt_result(public_key);
  ctxt_result = FHEMatMultMain(context, ctxt1, ctxt2, m, l, n, D0, D1, slot_b, l0, l1, LeftRot, UpRot);
  ctxt_result = Replicate(ctxt_result, D0 / m, m, UpRot);
  ctxt_result = FHEMatMultMain(context, ctxt_result, ctxt2, m, l, n, D0, D1, slot_b, l0, l1, LeftRot, UpRot);
  ctxt_result = Replicate(ctxt_result, D0 / m, m, UpRot);
  ctxt_result = FHEMatMultMain(context, ctxt_result, ctxt2, m, l, n, D0, D1, slot_b, l0, l1, LeftRot, UpRot);
  ctxt_result = Replicate(ctxt_result, D0 / m, m, UpRot);
  ctxt_result = FHEMatMultMain(context, ctxt_result, ctxt2, m, l, n, D0, D1, slot_b, l0, l1, LeftRot, UpRot);
  ctxt_result = Replicate(ctxt_result, D0 / m, m, UpRot);
  HELIB_NTIMER_STOP(Mult_2023b);
  helib::printNamedTimer(std::cout, "Mult_2023b");
  // Decrypt
  helib::Ptxt<helib::BGV> plaintext_result(context);
  std::cout << ctxt_result.capacity() << std::endl;
  secret_key.Decrypt(plaintext_result, ctxt_result);
  printf_ptxt(plaintext_result, D0, D1, slot_b, m, n, p);
}

void matrix_nosquare_2023(int m, int l, int n)
{
  long p = 202001;
  unsigned long m0 = 101 * 5 * 5 * 5;
  unsigned long r = 1;
  int l_0 = l, m_0 = m, n_0 = n;
  unsigned long bits = 120;
  unsigned long c = 2;
  const int slot_b = 100;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m0
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  std::cout << "  Generating key-switching matrices..." << std::endl;
  HELIB_NTIMER_START(KeyGen);
  // HeapProfilerStart("test"); //开始监测

  secret_key.GenSecKey(); //
  int down = 2;
  int right = 11;
  secret_key.GenKeySWmatrix(1, 2, 0, 0);
  secret_key.GenKeySWmatrix(1, 11, 0, 0);
  secret_key.setKeySwitchMap();

  const helib::PubKey &public_key = secret_key;
  vector<int> factors = {1, 2, 4, 5, 10, 20, 25, 50, 100};
  int D0, D1;
  context.printout();
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(m, l, 10, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M2 is equal to: " << std::endl;
  std::vector<std::vector<long>> M2 = randmat(l, n, 10, p, "M2.txt");
  std::cout << std::endl;
  cout << "Plaintext multiplication: M1*M2 mod " << p << " is equal to: " << endl;
  plain_text_multiplication(M1, M2, m, l, n, p);
  if (m <= l)
  {
    D0 = l;
    D1 = l >= n ? l : n;
    std::vector<long> LeftRot, UpRot;
    LeftRot.push_back(1);
    UpRot.push_back(1);
    int l0 = slot_b / D0;
    int l1 = slot_b / D1;
    for (int i = 1; i < D1; i++)
    {
      LeftRot.push_back(Modpower(right, l1 * (D1 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(right, l1 * (D1 - i), m0), 0, 0);
    }
    for (int i = 1; i < D0; i++)
    {
      UpRot.push_back(Modpower(down, l0 * (D0 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(down, l0 * (D0 - i), m0), 0, 0);
    }
    secret_key.setKeySwitchMap();
    helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
    // Encode
    Encode(ptxt1, ptxt2, M1, M2, slot_b, l0, l1);
    // Encrypt

    helib::Ctxt ctxt1(public_key), ctxt2(public_key);
    public_key.Encrypt(ctxt1, ptxt1);
    public_key.Encrypt(ctxt2, ptxt2);
    //
    HELIB_NTIMER_START(Mult_2023a);
    std::vector<helib::Ptxt<helib::BGV>> U = std::vector<helib::Ptxt<helib::BGV>>(l, helib::Ptxt<helib::BGV>(public_key));
    for (int i = 0; i < l; i++)
    {
      for (int j = 0; j < m; j++)
      {
        int h = j;
        int r = (j + i) % l;
        U[i][h * l0 * slot_b + l1 * r] = 1;
      }
    }
    std::vector<helib::Ctxt> ctxt1_temp = std::vector<helib::Ctxt>(l, helib::Ctxt(public_key));
    for (int i = 0; i < l; i++)
    {
      helib::Ctxt temp = ctxt1;
      temp.multByConstant(U[i]);
      ctxt1_temp[i] = Replicate(temp, D1, 1, LeftRot);
    }
    helib::Ctxt ctxt_result = ctxt2;
    ctxt_result.multiplyBy(ctxt1_temp[0]);
    for (int i = 1; i < l; i++)
    {
      helib::Ctxt temp = ctxt2;
      temp.smartAutomorph(UpRot[i]);
      temp.multiplyBy(ctxt1_temp[i]);
      ctxt_result.addCtxt(temp);
    }
    ctxt2.smartAutomorph(UpRot[1]);
    HELIB_NTIMER_STOP(Mult_2023a);
    helib::printNamedTimer(std::cout, "Mult_2023a");
    // Decrypt
    helib::Ptxt<helib::BGV> plaintext_result(context);
    secret_key.Decrypt(plaintext_result, ctxt_result);
    printf_ptxt(plaintext_result, D0, D1, slot_b, m_0, n_0, p);
  }
  else
  {
    D0 = m;
    D1 = l >= n ? l : n;
    std::vector<long> LeftRot, UpRot;
    LeftRot.push_back(1);
    UpRot.push_back(1);
    int l0 = slot_b / D0;
    int l1 = slot_b / D1;
    for (int i = 1; i < D1; i++)
    {
      LeftRot.push_back(Modpower(right, l1 * (D1 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(right, l1 * (D1 - i), m0), 0, 0);
    }
    for (int i = 1; i < D0; i++)
    {
      UpRot.push_back(Modpower(down, l0 * (D0 - i), m0));
      secret_key.GenKeySWmatrix(1, Modpower(down, l0 * (D0 - i), m0), 0, 0);
    }
    secret_key.setKeySwitchMap();
    helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
    // Encode
    Encode(ptxt1, ptxt2, M1, M2, slot_b, l0, l1);
    // Encrypt

    helib::Ctxt ctxt1(public_key), ctxt2(public_key);
    public_key.Encrypt(ctxt1, ptxt1);
    public_key.Encrypt(ctxt2, ptxt2);
    //
    HELIB_NTIMER_START(Mult_2023a);
    std::vector<helib::Ptxt<helib::BGV>> U = std::vector<helib::Ptxt<helib::BGV>>(l, helib::Ptxt<helib::BGV>(public_key));
    for (int i = 0; i < l; i++)
    {
      for (int j = 0; j < m; j++)
      {
        int h = j;
        int r = (j + i) % l;
        U[i][h * l0 * slot_b + l1 * r] = 1;
      }
    }
    std::vector<helib::Ctxt> ctxt1_temp = std::vector<helib::Ctxt>(l, helib::Ctxt(public_key));
    for (int i = 0; i < l; i++)
    {
      helib::Ctxt temp = ctxt1;
      temp.multByConstant(U[i]);
      ctxt1_temp[i] = Replicate(temp, D1, 1, LeftRot);
    }
    ctxt2 = Replicate(ctxt2, m / l, l, UpRot);
    helib::Ctxt ctxt_result = ctxt2;
    ctxt_result.multiplyBy(ctxt1_temp[0]);
    for (int i = 1; i < l; i++)
    {
      helib::Ctxt temp = ctxt2;
      temp.smartAutomorph(UpRot[i]);
      temp.multiplyBy(ctxt1_temp[i]);
      ctxt_result.addCtxt(temp);
    }
    ctxt2.smartAutomorph(UpRot[1]);
    HELIB_NTIMER_STOP(Mult_2023a);
    helib::printNamedTimer(std::cout, "Mult_2023a");
    // Decrypt
    helib::Ptxt<helib::BGV> plaintext_result(context);
    secret_key.Decrypt(plaintext_result, ctxt_result);
    printf_ptxt(plaintext_result, D0, D1, slot_b, m_0, n_0, p);
  }
}

void example_2023(int m, int l, int n)
{
  long p = 3221319251;
  unsigned long m0 = 101 * 5 * 5 * 5*3;
  unsigned long r = 1;
  int l_0 = l, m_0 = m, n_0 = n;
  unsigned long bits = 500;
  unsigned long c = 2;
  const int slot_b = 100;
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m0).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m0
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  std::cout << "  Generating key-switching matrices..." << std::endl;
  HELIB_NTIMER_START(KeyGen);
  // HeapProfilerStart("test"); //开始监测

  secret_key.GenSecKey(); //
  int down = 2;
  int right = 11;
  secret_key.GenKeySWmatrix(1, 2, 0, 0);
  secret_key.GenKeySWmatrix(1, 11, 0, 0);
  secret_key.setKeySwitchMap();

  const helib::PubKey &public_key = secret_key;
  vector<int> factors = {1, 2, 4, 5, 10, 20, 25, 50, 100};
  int D0, D1;
  context.printout();
  std::cout << "The matrix M1 is equal to: " << std::endl;
  std::vector<std::vector<long>> M1 = randmat(m, l, 10, p, "M1.txt");
  std::cout << std::endl;
  std::cout << "The matrix M2 is equal to: " << std::endl;
  std::vector<std::vector<long>> M2 = randmat(l, n, 10, p, "M2.txt");
  std::cout << std::endl;
  cout << "Plaintext multiplication: M1*M2*M2*M2 mod " << p << " is equal to: " << endl;
  std::vector<std::vector<long>> M;
  M = plain_text_multiplication_output(M1, M2, m, l, n, p);
  M = plain_text_multiplication_output(M, M2, m, l, n, p);
  M = plain_text_multiplication_output(M, M2, m, l, n, p);
  M = plain_text_multiplication_output(M, M2, m, l, n, p);
  printf_matrix(M, m, n, log10(p) + 2, p);
  D0 = l;
  D1 = l >= n ? l : n;
  std::vector<long> LeftRot, UpRot;
  LeftRot.push_back(1);
  UpRot.push_back(1);
  int l0 = slot_b / D0;
  int l1 = slot_b / D1;
  for (int i = 1; i < D1; i++)
  {
    LeftRot.push_back(Modpower(right, l1 * (D1 - i), m0));
    secret_key.GenKeySWmatrix(1, Modpower(right, l1 * (D1 - i), m0), 0, 0);
  }
  for (int i = 1; i < D0; i++)
  {
    UpRot.push_back(Modpower(down, l0 * (D0 - i), m0));
    secret_key.GenKeySWmatrix(1, Modpower(down, l0 * (D0 - i), m0), 0, 0);
  }
  secret_key.setKeySwitchMap();
  helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
  // Encode
  Encode(ptxt1, ptxt2, M1, M2, slot_b, l0, l1);
  // Encrypt

  helib::Ctxt ctxt1(public_key), ctxt2(public_key);
  public_key.Encrypt(ctxt1, ptxt1);
  public_key.Encrypt(ctxt2, ptxt2);
  //
  HELIB_NTIMER_START(Mult_2023a);
  std::vector<helib::Ptxt<helib::BGV>> U = std::vector<helib::Ptxt<helib::BGV>>(l, helib::Ptxt<helib::BGV>(public_key));
  for (int i = 0; i < l; i++)
  {
    for (int j = 0; j < m; j++)
    {
      int h = j;
      int r = (j + i) % l;
      U[i][h * l0 * slot_b + l1 * r] = 1;
    }
  }
  std::vector<helib::Ctxt> ctxt1_temp = std::vector<helib::Ctxt>(l, helib::Ctxt(public_key));
  for (int i = 0; i < l; i++)
  {
    helib::Ctxt temp = ctxt1;
    temp.multByConstant(U[i]);
    ctxt1_temp[i] = Replicate(temp, D1, 1, LeftRot);
  }
  helib::Ctxt ctxt_result = helib::Ctxt(public_key);
  for (int i = 0; i < l; i++)
  {
    helib::Ctxt temp = ctxt2;
    temp.smartAutomorph(UpRot[i]);
    temp.multiplyBy(ctxt1_temp[i]);
    ctxt_result.addCtxt(temp);
  }
  for (int L = 0; L < 3; L++)
  {
    for (int i = 0; i < l; i++)
    {
      helib::Ctxt temp = ctxt_result;
      temp.multByConstant(U[i]);
      ctxt1_temp[i] = Replicate(temp, D1, 1, LeftRot);
    }
    ctxt_result = helib::Ctxt(public_key);
    for (int i = 0; i < l; i++)
    {
      helib::Ctxt temp = ctxt2;
      temp.smartAutomorph(UpRot[i]);
      temp.multiplyBy(ctxt1_temp[i]);
      ctxt_result.addCtxt(temp);
    }
  }
  HELIB_NTIMER_STOP(Mult_2023a);
  helib::printNamedTimer(std::cout, "Mult_2023a");
  // Decrypt
  helib::Ptxt<helib::BGV> plaintext_result(context);
  secret_key.Decrypt(plaintext_result, ctxt_result);
  printf_ptxt(plaintext_result, D0, D1, slot_b, m, n, p);
}

int main()
{
  int m, l, n, thread;
  cout << "            ------------------------------- " << endl;
  cout << "             Momomorpic Matrix Operations" << endl;
  cout << "            ------------------------------- " << endl;
  cout << "Please input the size of matrix (m*l*n):" << endl;
  cin >> m;
  cin >> l;
  cin >> n;
  cout << "Please input the number of test:" << endl;
  int k;
  cin >> k;
  for (int i = 0; i < k; i++)
  {
    matrix_nosquare(m, l, n);
    // matrix_nosquare(20, 100, 5);
    matrix_nosquare_2023(m, l, n);
  }
  // example(5, 50, 50);
  // example_2023(5, 50, 50);
  return 0;
}
