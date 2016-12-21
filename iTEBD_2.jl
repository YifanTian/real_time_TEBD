function dosvdtrunc(AA,m)		# AA a matrix;  keep at most m states
    (u,d,v) = svd(AA)
    prob = dot(d,d)		# total probability
    mm = min(m,length(d))	# number of states to keep
    #=
    Dsum = 0
    sVn = 0
    n = 0
    #Dm = length(D)
      while (prob-Dsum)>1E-6 #&& n>=m
        n+=1
        Dsum += d[n]^2
        sVn += -d[n]^2*log(d[n]^2)
      end
    mm = max(n,md)
    =#
    d = d[1:mm]
    #d = d[1:mm]			# middle matrix in vector form
    trunc = prob - dot(d,d)
    sVn = sum(j->d[j]^2,1:mm)
    U = u[:,1:mm]
    V = v[:,1:mm]'
    (U,d,V,trunc,sVn)		# AA == U * diagm(d) * V	with error trunc
end

#=
function dosvdleftright(AA,m,toright)
    (U,d,V,trunc) = dosvdtrunc(AA,m)
    if toright
	V = diagm(d) * V
    else
	U = U * diagm(d)
    end
    (U,V,trunc)
end
=#

function dosvd4(AA,m,Dl)	# AA is ia * 2 * 2 * ib;  svd down the middle;  return two parts
    ia = size(AA,1)
    ib = size(AA,4)
    AA = reshape(AA,ia*2,2*ib)
    #(U,V,trunc) = dosvdleftright(AA,m,toright)
    (U,D,V,trunc,sVn) = dosvdtrunc(AA,m)
    mm = size(U,2)
    U = reshape(U,ia,2,mm)
    U = reshape(inv(Dl)*reshape(U,ia,2*mm),ia,2,mm)
    V = reshape(V,mm,2,ib)
    V = reshape(reshape(V,2*mm,ib)*inv(Dl),mm,2,ib)
    (U,diagm(D),V,trunc,sVn)
end

using TensorOperations
using PyPlot

sz = Float64[0.5 0; 0 -0.5]
sp = Float64[0 1; 0 0]
sm = sp'
Htwosite = Float64[sz[s1,s1p] * sz[s2,s2p] + 0.5 * (sp[s1,s1p] * sm[s2,s2p] + sm[s1,s1p] * sp[s2,s2p])
		    for s1=1:2, s2=1:2, s1p=1:2, s2p=1:2]
#Htwosite = Float64[sz[s1,s1p] * sz[s2,s2p] + 0.5 * (sp[s1,s1p] * sm[s2,s2p] + sm[s1,s1p] * sp[s2,s2p])
#		    for s1=1:2, s1p=1:2, s2=1:2, s2p=1:2]

tau = 0.1im
taugate = reshape(expm(-tau * reshape(Htwosite,4,4)),2,2,2,2)
#@show taugate

n = 2
#n = 28
#  Make initial product state in up down up down up down pattern (Neel state)
# Make first tensor a 1 x 2 x m tensor; and last is m x 2 x 1  (rather than vectors)
ps = [zeros(1,2,1) for i=1:n]
for i=1:n
    ps[i][1,iseven(i) ? 2 : 1,1] = 1.0
end
m = 10

A = ps[1]
B = ps[2]
D1 = zeros(1,1)
D1[1,1] = 1
D2 = zeros(1,1)
D2[1,1] = 1

nsweep = 1000
tstep = [1:1:nsweep]
SzAt = zeros(1,1000)
SzBt = zeros(1,1000)
trunct = zeros(1,1000)
vSnt = zeros(1,1000)

#for swp = 1:1000
#for swp = 1:1
    #@tensor begin
    @tensor	    AA[ml1,s1,s2,mr2] := D1[ml1,ml2] * A[ml2,s1,m1] * D2[m1,m2] * B[m2,s2,mr1] * D1[mr1,mr2]
    #@tensor	    AA[s1,s2,mr2,ml1] := D1[ml1,ml2] * A[ml2,s1,m1] * D2[m1,m2] * B[m2,s2,mr1] * D1[mr1,mr2]

    #@tensor      AA[s1p,s2p,mr2,ml1] := D1[ml1,ml2] * A[ml2,s1,m1] * D2[m1,m2] * B[m2,s2,mr1] * D1[mr1,mr2] * taugate[s1,s2,s1p,s2p]
    @tensor      AA[ml1,s1p,s2p,mr2] := D1[ml1,ml2] * A[ml2,s1,m1] * D2[m1,m2] * B[m2,s2,mr1] * D1[mr1,mr2] * taugate[s1,s2,s1p,s2p]
    @tensor 	    nor = scalar(AA[ml1,s1p,s2p,mr2] * AA[ml1,s1p,s2p,mr2])
    #end


      AA = reshape(AA,1*2,2*1)
      (U,D,V,trunc,sVn) = dosvdtrunc(AA,m)
      diagm(D)

    	AA *= 1.0 / sqrt(nor)
    	(A,D2,B,trunc,sVn) = dosvd4(AA,m,D1)

      @tensor      AA[s1p,s2p,mr2,ml1] := D1[ml1,ml2] * A[ml2,s1,m1] * D2[m1,m2] * B[m2,s2,mr1] * D1[mr1,mr2] * taugate[s1,s2,s1p,s2p]




      @tensor begin
    	    #AA[a,f,g,e] := Ai[a,b,c] * Ai1[c,d,e] * taugate[b,d,f,g]
    @tensor    AA[s2,s1,mr2,ml1] := D2[ml1,ml2] * B[ml2,s2,m1] * D1[m1,m2] * A[m2,s1,mr1] * D2[mr1,mr2]


      @tensor    AA[ml1,s2p,s1p,mr2] := D2[ml1,ml2] * B[ml2,s2,m1] * D1[m1,m2] * A[m2,s1,mr1] * D2[mr1,mr2] * taugate[s1,s2,s1p,s2p]

      @tensor    AA[s2p,s1p,ml1,mr2] := D2[ml1,ml2] * B[ml2,s2,m1] * D1[m1,m2] * A[m2,s1,mr1] * D2[mr1,mr2] * taugate[s1,s2,s1p,s2p]

    @tensor  nor = scalar(AA[ml1,s2p,s1p,mr2] * AA[ml1,s2p,s1p,mr2])
    	end
    	AA *= 1.0 / sqrt(nor)
    	(B,D1,A,trunc,vSn) = dosvd4(AA,m,D2)

      #dagA = conj!(A)
      dagA = conj!(A)
      @tensor SzA = scalar(D1[ml,mlrd]*A[mlrd,s1,m1d]*D2[m1d,m2]*D1[ml,mlru]*dagA[mlru,s1p,m1u]*D2[m1u,m2]*sz[s1,s1p])
      #dagB = conj!(B)
      dagB = conj(B)
      @tensor SzB = scalar(D2[ml,mlrd]*B[mlrd,s2,m1d]*D1[m1d,m2]*D2[ml,mlru]*dagB[mlru,s2p,m1u]*D1[m1u,m2]*sz[s2,s2p])

      #=
      #if swp%1000 == 0
      #  n += 1
        SzAt[swp] = real(SzA)
        SzBt[swp] = real(SzB)
        trunct[swp] = trunc
        vSnt[swp] = vSn
        @show SzA,SzB,trunc,sVn
      #end
      =#

      #@show size(A),size(B),trunc
      #@show SzA,SzB,trunc,sVn
#end

#plot(tstep[1:1000], SzAt[1:1000], color="red", linewidth=2.0, linestyle="--")
#plot(tstep[1:1000], SzBt[1:1000], color="red", linewidth=2.0, linestyle="--")

#plot(tstep[1:1000], trunct[1:1000], color="red", linewidth=2.0, linestyle="--")
#plot(tstep[1:1000], vSnt[1:1000], color="red", linewidth=2.0, linestyle="--")
