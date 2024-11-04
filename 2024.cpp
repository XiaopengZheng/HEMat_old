
#include "2024.h"
#include <helib/helib.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <helib/helib.h>
#include <vector>
#include <fstream>

using namespace std;
using namespace NTL;

helib::Ctxt Replicate(helib::Ctxt ctxt, int r, int s, std::vector<long> Rot)
{
  if (r == 1)
  {
    return ctxt;
  }
  if (r % 2 == 0)
  {
    helib::Ctxt temp = Replicate(ctxt, r / 2, s, Rot);
    helib::Ctxt temp2 = temp;
    temp2.smartAutomorph(Rot[(r / 2) * s]);
    temp.addCtxt(temp2);
    return temp;
  }
  else
  {
    helib::Ctxt temp = Replicate(ctxt, (r - 1) / 2, s, Rot);
    helib::Ctxt temp2 = temp;
    helib::Ctxt temp3 = ctxt;
    temp2.smartAutomorph(Rot[((r - 1) / 2) * s]);
    temp3.smartAutomorph(Rot[(r - 1) * s]);
    temp.addCtxt(temp2);
    temp.addCtxt(temp3);
    return temp;
  }
}

void printf_ptxt(helib::Ptxt<helib::BGV> &ptxt, int D0, int D1, int slot_b, int m, int n, long p)
{
  int l0 = slot_b / D0;
  int l1 = slot_b / D1;
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (long(ptxt[i * l0 * slot_b + l1 * j]) < p / 2)
        cout << setw(log10(p) + 2) << long(ptxt[i * l0 * slot_b + l1 * j]) << "";
      else
        cout << setw(log10(p) + 2) << long(ptxt[i * l0 * slot_b + l1 * j]) - p << "";
    }
    std::cout << std::endl;
  }
}

void printf_ptxt_block(std::vector<helib::Ptxt<helib::BGV>> &ptxt, int n, int slot, long p, int b_size)
{
  std::vector<long> M(n * n);
  const int l = 2;
  int b = n / b_size;
  for (int L = 0; L < b * b; L++)
  {
    int s = L / b;
    int t = L % b;
    for (int k = 0; k < 1; k++)
    {
      for (int i = 0; i < b_size; i++)
      {
        for (int j = 0; j < b_size; j++)
        {
          M[n * (b_size * s + i) + (b_size * t + j)] = long(ptxt[L][(b_size * b_size * k + i * b_size + j) * l]);
        }
      }
    }
  }

  std::ofstream fout;
  int w = int(log10(p)) + 2;
  fout.open("M1*M2_cypthertext.txt", ios::out);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      long a = long(M[i * n + j]) % p;
      if (a < p / 2)
      {
        fout << setw(w) << a << "";
      }
      else
      {
        fout << setw(w) << a - p << "";
      }
    }
    fout << endl;
  }
  cout << "See M1*M2_cypthertext.txt" << endl;
  fout.close();
  cout << endl;
}

int Modpower(int a, int b, int m)
{
  if (b == 1)
  {
    return a;
  }
  else
  {
    if (b != 1)
    {
      return (a * Modpower(a, b - 1, m)) % m;
    }
    else
    {
      int t = Modpower(a, b - 1, m);
      return (t * t) % m;
    }
  }
}

void plain_text_multiplication(vector<vector<long>> M1, vector<vector<long>> M2, int m, int l, int n, long p)
{
  vector<vector<long>> M(m, vector<long>(n));
  // long p = INT32_MAX;
  std::ofstream fout;
  fout.open("M1*M2_plaintext.txt", ios::out);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      M[i][j] = 0;
      for (int k = 0; k < l; k++)
      {
        M[i][j] = M[i][j] + M1[i][k] * M2[k][j];
      }
      long a = M[i][j] % p;
      if (abs(a) > p / 2)
      {
        a = (a > 0 ? a - p : a + p);
      }
      std::cout << setw(log10(p) + 2) << a << " ";
      fout << setw(log10(p) + 2) << a << " ";
    }
    std::cout << endl;
    fout << endl;
  }
  cout << "See M1*M2_plaintext.txt" << endl;
  fout.close();
  cout << endl;
}

vector<vector<long>> plain_text_multiplication_output(vector<vector<long>> M1, vector<vector<long>> M2, int m, int l, int n, long p)
{
  vector<vector<long>> M(m, vector<long>(n));
  // long p = INT32_MAX;
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      M[i][j] = 0;
      for (int k = 0; k < l; k++)
      {
        M[i][j] = M[i][j] + M1[i][k] * M2[k][j];
      }
    }
  }
  return M;
}

void printf_matrix(const std::vector<std::vector<long>> &M, int m, int n, int w, long &p)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      std::cout << setw(w) << M[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

void printf_matrix(const vector<vector<long>> &M, int m, int n, int w, long &p, ofstream &fout)
{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int a;
      if (M[i][j] % p < p / 2)
        a = M[i][j] % p;
      else
        a = M[i][j] % p - p;
      fout << setw(w) << a << " ";
    }
    fout << std::endl;
  }
}

std::vector<std::vector<long>> randmat(int m, int n, int q, long p, std::string s)
{
  std::vector<std::vector<long>> M(m, std::vector<long>(n, 0));
  int d = int(log10(q) + 2);
  srand((unsigned)time(NULL) + m);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      long a = rand() % q;
      if (a > q / 2)
      {
        a = a - q;
      }
      M[i][j] = a;
    }
  }
  if (n <= 64)
  {
    printf_matrix(M, m, n, d, p);
  }
  std::ofstream fout;
  fout.open(s, ios::out);
  printf_matrix(M, m, n, d, p, fout);
  cout << "See " << s << endl;
  fout.close();
  return M;
}

helib::Ctxt RotateAlign(helib::Context &context, helib::Ctxt ctxt, int l, int D0, int D1, int slot_b, int l0, int l1, bool t, std::vector<long> LeftRot, std::vector<long> UpRot)
{
  if (t == 1)
  {
    helib::Ctxt result = ctxt;
    std::vector<helib::Ptxt<helib::BGV>> U = std::vector<helib::Ptxt<helib::BGV>>(l, helib::Ptxt<helib::BGV>(context));
    for (int i = 0; i < l; i++)
    {
      for (int k = 0; k < D0; k++)
      {
        if (k % l == i)
        {
          for (int j = 0; j < D1; j++)
          {
            U[i][k * l0 * slot_b + l1 * j] = 1;
          }
        }
      }
    }
    result.multByConstant(U[0]);
    for (int i = 1; i < l; i++)
    {
      helib::Ctxt temp = ctxt;
      temp.multByConstant(U[i]);
      temp.smartAutomorph(LeftRot[i]);
      result.addCtxt(temp);
    }
    return result;
  }
  else
  {
    helib::Ctxt result = ctxt;
    std::vector<helib::Ptxt<helib::BGV>> V = std::vector<helib::Ptxt<helib::BGV>>(l, helib::Ptxt<helib::BGV>(context));
    for (int i = 0; i < l; i++)
    {
      for (int k = 0; k < D1; k++)
      {
        if (k % l == i)
        {
          for (int j = 0; j < D0; j++)
          {
            V[i][j * l0 * slot_b + l1 * k] = 1;
          }
        }
      }
    }
    result.multByConstant(V[0]);
    for (int i = 1; i < l; i++)
    {
      helib::Ctxt temp = ctxt;
      temp.multByConstant(V[i]);
      temp.smartAutomorph(UpRot[i]);
      result.addCtxt(temp);
    }
    return result;
  }
}

helib::Ctxt FHEMatMultMain(helib::Context &context, helib::Ctxt ctxt1, helib::Ctxt ctxt2, int m, int l, int n, int D0, int D1, int slot_b, int l0, int l1, std::vector<long> LeftRot, std::vector<long> UpRot)
{
  ctxt1 = RotateAlign(context, ctxt1, l, D0, D1, slot_b, l0, l1, 1, LeftRot, UpRot);
  ctxt2 = RotateAlign(context, ctxt2, l, D0, D1, slot_b, l0, l1, 0, LeftRot, UpRot);
  helib::Ctxt ctxt = ctxt1;
  ctxt.multiplyBy(ctxt2);
  int min;
  if (m<=l){
    min = n<m?n:m;
  }
  else{
    min = n<l?n:l;
  }
  for (int i = 1; i < min; i++)
  {
    helib::Ctxt temp1 = ctxt1;
    helib::Ctxt temp2 = ctxt2;
    temp1.smartAutomorph(LeftRot[i]);
    temp2.smartAutomorph(UpRot[i]);
    temp1.multLowLvl(temp2);
    ctxt.addCtxt(temp1);
  } 
  ctxt.reLinearize(); 
  return ctxt;
}

void Encode(helib::Ptxt<helib::BGV> &ptxt1, helib::Ptxt<helib::BGV> &ptxt2, std::vector<std::vector<long>> M1, std::vector<std::vector<long>> M2, int slot_b, int l0, int l1)
{
  int m = M1.size();
  int l = M1[0].size();
  int n = M2[0].size();
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < l; j++)
    {
      ptxt1[i * l0 * slot_b + l1 * j] = M1[i][j];
    }
  }

  for (int i = 0; i < l; i++)
  {
    for (int j = 0; j < n; j++)
    {
      ptxt2[i * l0 * slot_b + l1 * j] = M2[i][j];
    }
  }
}

void Encode(helib::Ptxt<helib::BGV> &ptxt1, std::vector<std::vector<long>> M1, int slot_b, int l0, int l1)
{
  int m = M1.size();
  int l = M1[0].size();
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < l; j++)
    {
      ptxt1[i * l0 * slot_b + l1 * j] = M1[i][j];
    }
  }
}

std::vector<helib::Ptxt<helib::BGV>> genU(helib::Context &context, int n, int baby, int giant, int slot)
{
  std::vector<helib::Ptxt<helib::BGV>> U = std::vector<helib::Ptxt<helib::BGV>>(2 * n - 1, helib::Ptxt<helib::BGV>(context));
  int s = 0;
  int m = slot / (n * n);
  for (int j = 0; j < baby - 1; j++)
  {
    for (int k = 0; k < j + 1; k++)
    {
      U[s][m * n * (n - j - 1) + m * (k + (n - 1))] = 1;
    }
    s++;
  }
  for (int i = -(giant - 1); i < 0; i++)
  {
    for (int j = 0; j < baby; j++)
    {
      int l = i * baby + j;
      for (int k = 0; k < n + l; k++)
      {
        U[s][m * n * (-l) + m * (k + ((-i) * baby))] = 1;
      }
      s++;
      //           std::cout << s - 1 << endl;
      // printf_ptxt(U[s - 1], 16, 1024 * 8, 1000);
      // std::cout << std::endl;
    }
  }

  for (int k = 0; k < n; k++)
  {
    U[s][m * k] = 1;
  }
  s++;

  // std::cout << s - 1 << endl;
  // printf_ptxt(U[s - 1], 16, 1024 * 8, 1000);
  // std::cout << std::endl;

  for (int i = 0; i < giant; i++)
  {
    for (int j = 0; j < baby; j++)
    {
      if (i != 0 || j != 0)
      {
        int l = i * baby + j;
        for (int k = n - 1; k >= l; k--)
        {
          U[s][m * n * (n - l) + m * (k - i * baby)] = 1;
        }
        s++;
        // std::cout << s - 1 << endl;
        // printf_ptxt(U[s - 1], 16, 1024 * 8, 1000);
        // std::cout << std::endl;
      }
    }
  }
  return U;
}

std::vector<helib::Ptxt<helib::BGV>> genV(helib::Context &context, int n, int slot)
{
  std::vector<helib::Ptxt<helib::BGV>> V;
  int m = slot / (n * n);
  for (int i = 0; i < n; i++)
  {
    helib::Ptxt<helib::BGV> ptxt(context);
    for (int j = 0; j < n; j++)
    {
      ptxt[m * n * j + m * i] = 1;
    }
    V.push_back(ptxt);
  }
  // for (int i =0;i<V.size();i++){
  //   printf_ptxt(V[i],n,slot,1000000);
  // }
  return V;
}

std::vector<std::vector<helib::Ptxt<helib::BGV>>> genW(helib::Context &context, int n, int slot)
{
  std::vector<std::vector<helib::Ptxt<helib::BGV>>> W;
  std::vector<helib::Ptxt<helib::BGV>> W1;
  std::vector<helib::Ptxt<helib::BGV>> W2;
  int m = slot / (n * n);
  for (int i = 0; i < n; i++)
  {
    helib::Ptxt<helib::BGV> ptxt1(context), ptxt2(context);
    for (int k = 0; k < n; k++)
    {
      for (int j = 0; j < i; j++)
      {
        ptxt1[m * n * k + m * j] = 1;
      }
      for (int j = i; j < n; j++)
      {
        ptxt2[m * n * k + m * j] = 1;
      }
    }
    W1.push_back(ptxt1);
    W2.push_back(ptxt2);
  }
  W.push_back(W1);
  W.push_back(W2);
  return W;
}

helib::Ctxt genA_BSGS(helib::Ctxt ctxt, std::vector<helib::Ptxt<helib::BGV>> U, int baby, int giant, int n, std::vector<long> LeftRotKey, std::vector<long> RightRotKey, const helib::PubKey &public_key)
{
  int m = baby * giant;
  helib::Ctxt ctxt_result(public_key);
  std::vector<helib::Ctxt> ctxt_baby = std::vector<helib::Ctxt>(baby, ctxt);
#pragma omp parallel for
  for (int j = 1; j < baby; j++)
  {
    ctxt_baby[j].smartAutomorph(RightRotKey[j - 1]);
    ctxt_baby[j].dropSmallAndSpecialPrimes();
  }
  std::vector<helib::Ctxt> ctxt_temp1 = ctxt_baby;
#pragma omp parallel for
  for (int j = 0; j < baby - 1; j++)
  {
    ctxt_temp1[j].multByConstant(U[j]);
  }

  for (int j = 0; j < baby - 1; j++)
  {
    ctxt_result.addCtxt(ctxt_temp1[j]);
  }
  ctxt_result.smartAutomorph(LeftRotKey[n - 2]);

  for (int i = 1; i < giant; i++)
  {
    helib::Ctxt ctxt_giant(public_key);
    std::vector<helib::Ctxt> ctxt_temp2 = std::vector<helib::Ctxt>(baby, ctxt);
#pragma omp parallel for
    for (int j = 0; j < baby; j++)
    {
      ctxt_temp2[j] = ctxt_baby[j];
      ctxt_temp2[j].multByConstant(U[j + i * baby - 1]);
    }
    for (int j = 0; j < baby; j++)
      ctxt_giant.addCtxt(ctxt_temp2[j]);
    ctxt_giant.smartAutomorph(LeftRotKey[(giant - i) * baby - 1]);
    ctxt_result.addCtxt(ctxt_giant);
  }

  for (int i = 0; i < giant; i++)
  {
    helib::Ctxt ctxt_giant(public_key);
    std::vector<helib::Ctxt> ctxt_temp2 = std::vector<helib::Ctxt>(baby, ctxt);
#pragma omp parallel for
    for (int j = 0; j < baby; j++)
    {
      ctxt_temp2[j] = ctxt_baby[j];
      ctxt_temp2[j].multByConstant(U[i * baby + j + n - 1]);
      // std::cout << RightRotKey[i - 1] << " " << i + n - 1 << " "
    }
    for (int j = 0; j < baby; j++)
    {
      ctxt_giant.addCtxt(ctxt_temp2[j]);
    }
    if (i > 0)
    {
      ctxt_giant.smartAutomorph(RightRotKey[i * baby - 1]);
    }
    ctxt_result.addCtxt(ctxt_giant);
  }
  return ctxt_result;
}

helib::Ctxt genB_BSGS(helib::Ctxt ctxt, std::vector<helib::Ptxt<helib::BGV>> V_BSGS, int baby, int giant, std::vector<long> BlockRotKey, const helib::PubKey &public_key)
{
  helib::Ctxt ctxt_result = ctxt;
  std::vector<helib::Ctxt> ctxt_baby = std::vector<helib::Ctxt>(baby, ctxt);
  ctxt_baby[0] = ctxt;
#pragma omp parallel for
  for (int j = 1; j < baby; j++)
  {
    helib::Ctxt ctxt_temp = ctxt;
    ctxt_temp.smartAutomorph(BlockRotKey[j - 1]);
    ctxt_baby[j] = ctxt_temp;
  }
  for (int i = 0; i < giant; i++)
  {
    helib::Ctxt ctxt_giant(public_key);
    std::vector<helib::Ctxt> ctxt_temp2 = std::vector<helib::Ctxt>(baby, ctxt);
#pragma omp parallel for
    for (int j = 0; j < baby; j++)
    {
      ctxt_temp2[j] = ctxt_baby[j];
      ctxt_temp2[j].multByConstant(V_BSGS[baby * i + j]);
    }
    for (int j = 0; j < baby; j++)
    {
      ctxt_giant.addCtxt(ctxt_temp2[j]);
    }
    if (i == 0)
    {
      ctxt_result = ctxt_giant;
    }
    else
    {
      ctxt_giant.smartAutomorph(BlockRotKey[baby * i - 1]);
      ctxt_result.addCtxt(ctxt_giant);
    }
  }
  return ctxt_result;
}

helib::Ctxt multStep3(helib::Ctxt ctxt1, helib::Ctxt ctxt2, int n, std::vector<long> LeftRotKey, std::vector<long> RightRotKey, std::vector<long> BlockRotKey, std::vector<std::vector<helib::Ptxt<helib::BGV>>> W)
{
  ctxt1.dropSmallAndSpecialPrimes();
  ctxt2.dropSmallAndSpecialPrimes();
  helib::Ctxt ctxt_result = ctxt1;
  std::vector<helib::Ctxt> ctxt_temp2 = std::vector<helib::Ctxt>(n, ctxt2);
#pragma omp parallel for
  for (int i = 0; i < n; i++)
  {
    if (i == 0)
    {
      ctxt_result.multiplyBy(ctxt2);
    }
    else
    {
      helib::Ctxt ctxt_temp11 = ctxt1;
      helib::Ctxt ctxt_temp12 = ctxt1;
      ctxt_temp11.multByConstant(W[0][i]);
      ctxt_temp12.multByConstant(W[1][i]);
      ctxt_temp11.smartAutomorph(RightRotKey[n - i - 1]);
      ctxt_temp12.smartAutomorph(LeftRotKey[i - 1]);
      ctxt_temp11.addCtxt(ctxt_temp12);
      ctxt_temp2[i].smartAutomorph(BlockRotKey[i - 1]);
      ctxt_temp2[i].multLowLvl(ctxt_temp11);
    }
  }
  for (int i = 1; i < n; i++)
  {
    ctxt_result.addCtxt(ctxt_temp2[i]);
  }
  ctxt_result.reLinearize();
  return ctxt_result;
}

helib::Ctxt MatrixMult(helib::Ctxt ctxt1, helib::Ctxt ctxt2, std::vector<helib::Ptxt<helib::BGV>> U, std::vector<helib::Ptxt<helib::BGV>> V, std::vector<std::vector<helib::Ptxt<helib::BGV>>> W, std::vector<long> LeftRotKey, std::vector<long> RightRotKey, std::vector<long> BlockRotKey, int n, int baby, int giant, const helib::PubKey &public_key)
{
  helib::Ctxt ctxt_result = ctxt1;
  // HELIB_NTIMER_START(Total);
  //  Generating B
  // HELIB_NTIMER_START(Gen_A);
  ctxt1 = genA_BSGS(ctxt1, U, baby, giant, n, LeftRotKey, RightRotKey, public_key);
  // HELIB_NTIMER_STOP(Gen_A);
  helib::printNamedTimer(std::cout, "Gen_A");

  // HELIB_NTIMER_START(Gen_B);
  ctxt2 = genB_BSGS(ctxt2, V, baby, giant, BlockRotKey, public_key);
  // HELIB_NTIMER_STOP(Gen_B);
  // helib::printNamedTimer(std::cout, "Gen_B");

  // HELIB_NTIMER_START(Step3);
  ctxt_result = multStep3(ctxt1, ctxt2, n, LeftRotKey, RightRotKey, BlockRotKey, W);
  // HELIB_NTIMER_STOP(Step3);
  // helib::printNamedTimer(std::cout, "Step3");

  // HELIB_NTIMER_STOP(Total);
  // helib::printNamedTimer(std::cout, "Total");
  return ctxt_result;
}
