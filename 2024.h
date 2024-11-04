#include <helib/helib.h>
#include <vector>
#include <omp.h>
#include <cmath>
#include <gperftools/heap-profiler.h>
#include <omp.h>
#include <stdlib.h>
#include <string>

struct Ctxt_leval
{
   helib::Ctxt m_ctxt;
   int m_L;
   Ctxt_leval(helib::Ctxt ctxt, int L): m_ctxt(ctxt), m_L(L){};
};


void Encode(helib::Ptxt<helib::BGV> &ptxt1, helib::Ptxt<helib::BGV> & ptxt2,std::vector<std::vector<long>>  M1, std::vector<std::vector<long>>  M2,int slot_b,int l0,int l1);

void Encode(helib::Ptxt<helib::BGV> &ptxt1, std::vector<std::vector<long>>  M1, int slot_b,int l0,int l1);

void printf_ptxt(helib::Ptxt<helib::BGV> &ptxt, int D0, int D1, int slot_b, int m, int n, long p);

void printf_ptxt_block(std::vector<helib::Ptxt<helib::BGV>> &ptxt, int n, int slot, long p, int b_size);

int Modpower(int a,int b,int m);

void plain_text_multiplication(std::vector<std::vector<long>> M1, std::vector<std::vector<long>> M2, int m, int l, int n, long p);

std::vector<std::vector<long>> plain_text_multiplication_output(std::vector<std::vector<long>> M1, std::vector<std::vector<long>> M2, int m, int l, int n, long p);

helib::Ctxt RotateAlign(helib::Context &context, helib::Ctxt ctxt, int l, int D0, int D1, int slot_b, int l0, int l1, bool t, std::vector<long> LeftRot, std::vector<long> UpRot);

helib::Ctxt FHEMatMultMain(helib::Context &context, helib::Ctxt ctxt1, helib::Ctxt ctxt2, int m, int l, int n, int D0, int D1, int slot_b, int l0, int l1, std::vector<long> LeftRot, std::vector<long> UpRot);

helib::Ctxt Replicate(helib::Ctxt ctxt, int r, int s,std::vector<long> Rot);

std::vector<std::vector<long>> randmat(int m, int n, int q, long p, std::string s);

void printf_matrix(const std::vector<std::vector<long>> &M, int m, int n, int w, long &p);

void printf_matrix(const std::vector<std::vector<long>> &M, int m, int n, int w, long &p, std::ofstream &fout);

std::vector<helib::Ptxt<helib::BGV>> genU(helib::Context &context, int n, int baby, int giant, int slot);

std::vector<helib::Ptxt<helib::BGV>> genV(helib::Context& context, int n,int slot);

std::vector<std::vector<helib::Ptxt<helib::BGV>>> genW(helib::Context& context, int n,int slot);

helib::Ctxt genA_BSGS(helib::Ctxt ctxt, std::vector<helib::Ptxt<helib::BGV>> U_BSGS, int baby, int giant, int n, std::vector<long> LeftRotKey, std::vector<long> RightRotKey,const helib::PubKey &public_key);

helib::Ctxt genB_BSGS(helib::Ctxt ctxt, std::vector<helib::Ptxt<helib::BGV>> V_BSGS, int baby, int giant, std::vector<long> BlockRotKey, const helib::PubKey &public_key);

helib::Ctxt multStep3(helib::Ctxt ctxt1, helib::Ctxt ctxt2, int n, std::vector<long> LeftRotKey, std::vector<long> RightRotKey, std::vector<long> BlockRotKey,std::vector<std::vector<helib::Ptxt<helib::BGV>>> W);

helib::Ctxt MatrixMult(helib::Ctxt ctxt1, helib::Ctxt ctxt2, std::vector<helib::Ptxt<helib::BGV>> U, std::vector<helib::Ptxt<helib::BGV>> V, std::vector<std::vector<helib::Ptxt<helib::BGV>>> W, std::vector<long> LeftRotKey, std::vector<long> RightRotKey, std::vector<long> BlockRotKey, int n, int baby, int giant, const helib::PubKey &public_key);


