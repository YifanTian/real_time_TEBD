#include "itensor/all.h"

using namespace itensor;
using std::vector;
//std::cin>>response;

char response;

vector<ITensor>
makeH(SiteSet const& sites)
    {
    auto N = sites.N();
    auto H = vector<ITensor>(N+1);
    for(auto b : range1(N-1))
        {
        //Make S.S operator on sites b, b+1
        H.at(b) = sites.op("Sz",b)*sites.op("Sz",b+1)
                   + 0.5*sites.op("S+",b)*sites.op("S-",b+1)
                   + 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        //Print(H.at(b));
        //PAUSE
        }
    return H;
    }

vector<ITensor>
makeGates(vector<ITensor> const& H,
          Real tau)
    {
    auto gates = H;
    for(auto& g : gates)
        {
        if(!g) continue;
        g = expHermitian(-tau*g);
        }
    return gates;
    }

void
doSVD(MPS & psi,
      int b,
      ITensor phi,
      Direction dir,
      Real cutoff)
    {
    auto U = psi.A(b);
    ITensor D,V;
    svd(phi,U,D,V,{"Cutoff",cutoff});
    if(dir == Fromleft)
        {
        //multiply D into V
        psi.setA(b,U);
        psi.setA(b+1,D*V);
        }
    else
        {
        //multiply D into U
        psi.setA(b,U*D);
        psi.setA(b+1,V);
        }
    }

ITensor
applyGate(ITensor phi,
          ITensor gate)
    {
    //TODO:

    phi = phi*gate;     //1. Apply gate to phi using * operator

    //println("phi:",phi);
    //std::cin>>response;

    phi.mapprime(1,0);      //2. Restore original prime level of phi's indices
    phi /= norm(phi);      //3. Normalize phi by dividing by norm(phi)

    return phi;
    }



int
main()
    {
    char response;

    int N = 4;
    Real tau = 0.001;
    Real cutoff = 1E-12;
    int nsweep = 2;
    auto Energybond = 0.0;

    auto sites = SpinHalf(N);

    auto state = InitState(sites);
    for(auto n : range1(N))
        {
        if(n%2==1) state.set(n,"Up");
        else       state.set(n,"Dn");
        }

    auto psi = MPS(state);

    auto H = makeH(sites);
    auto gates = makeGates(H,tau);

    auto U1 = psi.A(1);
    auto phi = psi.A(1)*psi.A(1+1);
    ITensor D1,D2,V1,U2;
    svd(phi,U1,D1,V1,{"Cutoff",cutoff});
    //println("here1");
    phi = psi.A(3)*psi.A(3+1);
    auto  V2 = psi.A(4);
    svd(phi,U2,D2,V2,{"Cutoff",cutoff});
    //println("here2");
    //PrintData(D1*psi.A(2));
    //PAUSE
    ITensor U,D,V;
    for(auto sw : range1(nsweep))
      {

        //ITensor U,D,V;
        //auto Energy = 0.0;
        //for(auto b : range1(N-1))   {
        //psi.position(b);
          println("here1");
        if(sw == 1) {
          phi = D1*psi.A(2)*psi.A(2+1)*D2; }
        else {
        phi = psi.A(2)*D*psi.A(2+1);
        //PrintData(psi.A(2)*D*psi.A(2+1));
        }
        //auto bra = prime(phi,Site);
        //Energy += (bra*H.at(b)*phi).real(); //measurement of energy <S.S>
        phi = applyGate(phi,gates.at(2));

        U = D1*psi.A(2);   //??
        svd(phi,U,D,V,{"Cutoff",cutoff});
        //doSVD(psi,b,phi,Fromleft,cutoff); }

        //D1 = D; D2 = D; //??                    //translational symmetry
        //D1inv = Dinversion(D1);
        //D2inv = Dinversion(D2);
        ITensor D1inv,D2inv;
        D1inv = D1;D2inv = D2;                  //translational symmetry
        //PrintData(D1*psi.A(2)*D*psi.A(2+1)*D2);
        PAUSE
        psi.setA(2,D1inv*U);
        psi.setA(2+1,V*D2inv);
        PrintData(psi.A(2)*D*psi.A(2+1));
        println("sw1");
        //D1 = D; D2 = D;

        //printfln("Half 1: Energy = %.10f, Energy per site = %.10f",Energy,Energy/N);

      }

        //if(b == N/2) {
          println("Bond dimension at center = ",commonIndex(psi.A(2),psi.A(2+1)).m());
        //}


      //  }

    return 0;
    }
