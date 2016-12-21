


using TensorOperations
using PyPlot


sz = Float64[0.5 0;0 -0.5]; sp = Float64[0 1;0 0]; sm = sp'
Htwosite = Float64[sz[s1,s1p]*sz[s2,s2p]+0.5*(sp[s1,s1p]*sm[s2,s2p]+sm[s1,s1p]*sp[s2,s2p])
          for s1 = 1:2, s2 = 1:2, s1p = 1:2, s2p = 1:2]
#tau = 0.01
tau = 0.01im
taugate = reshape(expm(-tau*reshape(Htwosite,4,4)),2,2,2,2) # s1,s1p,s2,s2p?
cutoff = 1E-6


A = zeros(1,2,1)
B = zeros(1,2,1)

A[1,1,1] = 1.0
B[1,2,1] = 1.0

@tensor SzA = scalar(A[ml,s1,mr]*sz[s1,s1p]*A[ml,s1p,mr])
@tensor SzB = scalar(B[ml,s1,mr]*sz[s1,s1p]*B[ml,s1p,mr])


D1 = [1]
D2 = [1]
D1 =reshape(D1,1,1)
D2 =reshape(D2,1,1)

@tensor SzA = scalar(D1[ml,mlrd]*A[mlrd,s1,m1d]*D2[m1d,m2]*D1[ml,mlru]*A[mlru,s1p,m1u]*D2[m1u,m2]*sz[s1,s1p])
@tensor SzB = scalar(D2[ml,mlrd]*B[mlrd,s2,m1d]*D1[m1d,m2]*D2[ml,mlru]*B[mlru,s2p,m1u]*D1[m1u,m2]*sz[s2,s2p])


#function twositephi(D1,A,D2,B)

function truncationUDV(AA)
  (a,b,c,d) = size(AA)
  AA = reshape(AA,a*b,c*d)
  (U,D,V) = svd(AA)
    Dsum = 0
    Dm = length(D)
    for i = 1:length(D)
      Dsum += D[i]^2
      if i>10
        if (1-Dsum)<1E-12
          Dm = i
          break
       end
      end
    end
  D = diagm(D[1:Dm])
  A = reshape(U[:,1:Dm],a,b*Dm)
  B = reshape(transpose(V[:,1:Dm]),Dm*c,d)
  #B = reshape(V,c*d*c,d)
  return (A,D,B,Dm)
end

#for sw = 1:nsweep
for sw = 1:1000
  @tensor phi1[ml,s1,s2,mr] := D1[ml,mlr]*A[mlr,s1,m1]*D2[m1,m2]*B[m2,s2,mrl]*D1[mrl,mr]
  @tensor AA[ml,s1p,s2p,mr] := phi1[ml,s1,s2,mr]*taugate[s1,s2,s1p,s2p]

  (a,b,c,d) = size(AA)
  (A,D2,B,Dm) = truncationUDV(AA)
  (ml1,ml2) = size(D1)
  A = reshape(inv(D1)*A,ml1,b,Dm)
  (mr1,mr2) = size(D1)
  B = reshape(B*inv(D1),Dm,c,mr2)

  @tensor SzA = scalar(D1[ml,mlrd]*A[mlrd,s1,m1d]*D2[m1d,m2]*D1[ml,mlru]*A[mlru,s1p,m1u]*D2[m1u,m2]*sz[s1,s1p])
  @tensor SzB = scalar(D2[ml,mlrd]*B[mlrd,s2,m1d]*D1[m1d,m2]*D2[ml,mlru]*B[mlru,s2p,m1u]*D1[m1u,m2]*sz[s2,s2p])
  @show size(A),size(B),size(D1),size(D2)

  @tensor phi2[ml,s1,s2,mr] := D2[ml,mlr]*B[mlr,s1,m1]*D1[m1,m2]*A[m2,s2,mrl]*D2[mrl,mr]
  @tensor AA[ml,s1p,s2p,mr] := phi2[ml,s1,s2,mr]*taugate[s1,s2,s1p,s2p]

  (a,b,c,d) = size(AA)
  (B,D1,A,Dm) = truncationUDV(AA)
  (ml1,ml2) = size(D2)
  B = reshape(inv(D2)*B,ml1,b,Dm)
  (mr1,mr2) = size(D2)
  A = reshape(A*inv(D2),Dm,c,mr2)

  @tensor SzA = scalar(D1[ml,mlrd]*A[mlrd,s1,m1d]*D2[m1d,m2]*D1[ml,mlru]*A[mlru,s1p,m1u]*D2[m1u,m2]*sz[s1,s1p])
  @tensor SzB = scalar(D2[ml,mlrd]*B[mlrd,s2,m1d]*D1[m1d,m2]*D2[ml,mlru]*B[mlru,s2p,m1u]*D1[m1u,m2]*sz[s2,s2p])
  @show SzA,SzB,(ml1,ml2),Dm
  #@show size(A),size(B),size(D1),size(D2)
end
#for sw = 1:nsweep

#end
