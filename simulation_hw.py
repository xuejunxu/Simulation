import math
# # first generator

import csv
with open('cimp_rv.csv', 'r') as f:
    reader = csv.reader(f)
    #print(reader)
    rv = list(reader)
# # chi-square test

print(len(rv))

num_rv=[]


n_rep=32768

k=16

z_0 = 0.84611665575369589265

d=3

lamda=1-10**(-11)


for i in range(n_rep):

	num_rv.append(int(str(rv[i+1])))

generating random variables, storing them in rv_storer list

rv_storer=[]

for i in range(n_rep):
	z_i=2*lamda*(0.5-abs(z_0-0.5))
	rv_storer.append(z_i)
	z_0=z_i


### starting chi-square test

# divide (0,1) to k subintervals

# create a k list and initialize all k intervals to 0

k_list=[0]*k

#print(k_list)

# use k_list to store Nj values

for j in range(k):
	for num in range(n_rep):
		if (rv_storer[num]>=(j-1)/k) and (rv_storer[num]<(j/k)):
			k_list[j]=k_list[j]+1


#print(len(k_list))

# calculate chi-square

chi_stat=0

for count in range(k):
	chi_stat=chi_stat+(k_list[count]-n_rep/k)**2

chi_stat=k/n_rep*chi_stat

print(chi_stat)

# if use alpha=0.05, calculate chi-square k-1, 1-alpha, z(1-alpha)=1.96



chi_approx=(k-1)*(1-2/(9*(k-1))-1.96*(2/(9*(k-1)))**0.5)**3

print(chi_approx)


z_alpha=(2/(9*(k-1)))**(-0.5)*(((1/(k-1))*(chi_stat))**(1/3)-1+2/(9*(k-1)))

print(z_alpha)


### ending chi-square test


### starting K-S test

# order statistics

rv_for_ks=rv_storer

rv_for_ks.sort()

print(rv_for_ks)

# dn_minus

dn_minus=[]
dn_plus=[]

for i in range(len(rv_for_ks)):
	#print(i)
	#print(i/n_rep-rv_for_ks[i])
	dn_minus.append(i/n_rep-rv_for_ks[i])
	dn_plus.append(rv_for_ks[i]-(i-1)/n_rep)

print(dn_minus)
print(dn_plus)

dn_p_max=max(dn_plus)
dn_m_max=max(dn_minus)

# print (dn_p_max)
# print (dn_m_max)

if (dn_p_max>=dn_m_max):
	dn=dn_p_max
else:
	dn=dn_m_max

print(dn)

ks_stat=(n_rep**0.5+0.12+0.11/(n_rep)**0.5)*dn_m_max

print(ks_stat)


### ending K-S test


### starting Serial test

# initial 2d list of U

serial_u = [[0 for i in range(d)] for j in range(n_rep)]

for i in range(n_rep):
	serial_u[i][0]=rv_storer[d*i]
	serial_u[i][1]=rv_storer[d*i+1]
	serial_u[i][2]=rv_storer[d*i+2]

print(serial_u)

# use 2d_list to store Nj values

# for j in range(k):
# 	for num in range(n_rep):
# 		if (rv_storer[num]>=(j-1)/k) and (rv_storer[num]<(j/k)):
# 			k_list[j]=k_list[j]+1

# # create a k list and initialize all k intervals to 0

N_j_d=[[[0 for i in range(k)] for j in range(k)] for m in range(k)]

print(len(N_j_d))

for j0 in range(k):
	for j1 in range(k):
		for j2 in range(k):
			for i in range(n_rep):
				if serial_u[i][0]>=(j0-1)/k and serial_u[i][0]<(j0/k) and serial_u[i][1]>=(j1-1)/k and serial_u[i][1]<(j1/k) and serial_u[i][2]>=(j2-1)/k and serial_u[i][2]<(j2/k):
					N_j_d[j0][j1][j2]+=1

print(N_j_d)

# calculate chi-stat

chi_d=0

for i in range(k):
	for j in range(k):
		for m in range(k):
			chi_d+=(N_j_d[i][j][m]-n_rep/(k**d))

chi_d=k**d/n_rep

print(chi_d)

## ending serial test

## starting runs test

a_matrix=[[4529.4,9044.9,13568,18091,22615,27892],[9044.9,18097,27139,36187,45234,55789],[13568,27139,40721,54281,67852,83685],[18091,36187,54281,72414,90470,111580],[22615,45234,67852,90470,113262,139476],[27892,55789,83685,111580,139476,172860]]

b_mat=[1/6,5/24,11/120,19/720,29/540,1/840]

r=[n_rep,0,0,0,0,0]

Rn=0

for i in range(6):
	for j in range(6):
		sum=a_matrix[i][j]*(r[i]-n_rep*b_mat[i])*(r[j]-n_rep*b_mat[j])
		Rn+=sum

Rn=1/n_rep*Rn

print(Rn)

## Ending runs test

## starting autocorrelation test

lag=6

h=math.floor(((n_rep-1)/lag)-1)

p_j=0

for i in range(h):
	sum_auto=rv_storer[i*lag]*rv_storer[(i+1)*lag]
	p_j+=sum_auto

p_j=12/(1+h)*p_j-3

var_pj=(13*h+7)/((h+1)**2)

Aj=p_j/math.sqrt(var_pj)

print(Aj)



## ending auto test

print(rv)
