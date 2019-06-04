import math
import turtle

### LCG
z_0 = 12345
a=1140671485
c=12820163
m=2**24

# n is the number of iterations

n=1000


lcg_rv1=[]
for i in range(n):
    z_i=(a*z_0+c)%m
    lcg_rv1.append(z_i/m)
    #print(z_i/m)
    z_0=z_i


z_0=23456

lcg_rv2=[]
for i in range(n):
    z_i=(a*z_0+c)%m
    lcg_rv2.append(z_i/m)
    #print(z_i/m)
    z_0=z_i

z_0=34567

lcg_rv3=[]
for i in range(n):
    z_i=(a*z_0+c)%m
    lcg_rv3.append(z_i/m)
    #print(z_i/m)
    z_0=z_i


### MRG 
n=1000

z1_0=1111
z1_1=2222
z1_2=3333

z2_0=4444
z2_1=5555
z2_2=6666

a1_1=0
a1_2=1403580
a1_3=-810728

a2_1=527612
a2_2=0
a2_3=-1370589

z1_list=[z1_0,z1_1,z1_2]
z2_list=[z2_0,z2_1,z2_2]

for i in range(3,n):
    z1_i=(a1_1*z1_list[i-1]+a1_2*z1_list[i-2]+a1_3*z1_list[i-3])%(2**32-209)
    z2_i=(a2_1*z1_list[i-1]+a2_2*z1_list[i-2]+a2_3*z1_list[i-3])%(2**32-22853)
    z1_list.append(z1_i)
    z2_list.append(z2_i)

mcg_rv1=[]
for i in range(n):
    y_i=(z1_list[i]-z2_list[i])%(2**32-209)
    if y_i>0:
        u_i=y_i/(2**32-208)
    else:
        u_i=(2**32-209)/(2**32-208)
    mcg_rv1.append(u_i)

## mcg 2

z2_2=1234

for i in range(3,n):
    z1_i=(a1_1*z1_list[i-1]+a1_2*z1_list[i-2]+a1_3*z1_list[i-3])%(2**32-209)
    z2_i=(a2_1*z1_list[i-1]+a2_2*z1_list[i-2]+a2_3*z1_list[i-3])%(2**32-22853)
    z1_list.append(z1_i)
    z2_list.append(z2_i)

mcg_rv2=[]
for i in range(n):
    y_i=(z1_list[i]-z2_list[i])%(2**32-209)
    if y_i>0:
        u_i=y_i/(2**32-208)
    else:
        u_i=(2**32-209)/(2**32-208)
    mcg_rv2.append(u_i)

## mcg 3

z2_1=3457

d=3

for i in range(3,n*d):
    z1_i=(a1_1*z1_list[i-1]+a1_2*z1_list[i-2]+a1_3*z1_list[i-3])%(2**32-209)
    z2_i=(a2_1*z1_list[i-1]+a2_2*z1_list[i-2]+a2_3*z1_list[i-3])%(2**32-22853)
    z1_list.append(z1_i)
    z2_list.append(z2_i)

mcg_rv3=[]
for i in range(n):
    y_i=(z1_list[i]-z2_list[i])%(2**32-209)
    if y_i>0:
        u_i=y_i/(2**32-208)
    else:
        u_i=(2**32-209)/(2**32-208)
    mcg_rv3.append(u_i)


# serial test mcg

for i in range(3,n*d):
    z1_i=(a1_1*z1_list[i-1]+a1_2*z1_list[i-2]+a1_3*z1_list[i-3])%(2**32-209)
    z2_i=(a2_1*z1_list[i-1]+a2_2*z1_list[i-2]+a2_3*z1_list[i-3])%(2**32-22853)
    z1_list.append(z1_i)
    z2_list.append(z2_i)

mcg_rv1_serial=[]
for i in range(n*d):
    y_i=(z1_list[i]-z2_list[i])%(2**32-209)
    #print("y_i is ",y_i)
    if y_i!=0:
        u_i=y_i/(2**32-208)
    else:
        u_i=(2**32-209)/(2**32-208)
    mcg_rv1_serial.append(u_i)

# sample mean

sum=0 

for i in range(n):
    sum=lcg_rv1[i]+sum
sample_mean_lcg1=sum/n
print("sample mean is", sample_mean_lcg1)

# sample std

sample_std_lcg1=0

for i in range(n):
    sample_std_lcg1=sample_std_lcg1+(lcg_rv1[i]-sample_mean_lcg1)**2
sample_std_lcg1=sample_std_lcg1/(n-1)
sample_std_lcg1=sample_mean_lcg1**0.5
print("sample standard deviation is",sample_std_lcg1)

# frequency chi-square test

k=100

### starting chi-square test

# divide (0,1) to k subintervals

# create a k list and initialize all k intervals to 0

k_list=[0]*k

#print(k_list)

# use k_list to store Nj values

# for j in range(k):
# 	for num in range(n):
# 		if (lcg_rv1[num]>=(j-1)/k) and (lcg_rv1[num]<(j/k)):
# 			k_list[j]=k_list[j]+1

# for MCG
for j in range(k):
    for num in range(n):
        if (mcg_rv1[num]>=(j-1)/k) and (mcg_rv1[num]<(j/k)):
            k_list[j]=k_list[j]+1

#print(len(k_list))

# calculate chi-square

chi_stat=0

for count in range(k):
	chi_stat=chi_stat+(k_list[count]-n/k)**2

chi_stat=k/n*chi_stat

print("chi statistic",chi_stat)

# if use alpha=0.05, calculate chi-square k-1, 1-alpha, z(1-alpha)=1.96

chi_approx=(k-1)*(1-2/(9*(k-1))-1.96*(2/(9*(k-1)))**0.5)**3

print("chi statistics approximation",chi_approx)


z_alpha=(2/(9*(k-1)))**(-0.5)*(((1/(k-1))*(chi_stat))**(1/3)-1+2/(9*(k-1)))

print("z score approximation",z_alpha)


### ending chi-square test

### starting the K-S test

# order statistics

#rv_for_ks=lcg_rv1

rv_for_ks=mcg_rv1

rv_for_ks.sort()

#print(rv_for_ks)

# dn_minus

dn_minus=[]
dn_plus=[]

for i in range(len(rv_for_ks)):
	#print(i)
	#print(i/n_rep-rv_for_ks[i])
	dn_minus.append(i/n-rv_for_ks[i])
	dn_plus.append(rv_for_ks[i]-(i-1)/n)

# print(dn_minus)
# print(dn_plus)

dn_p_max=max(dn_plus)
dn_m_max=max(dn_minus)

# print (dn_p_max)
# print (dn_m_max)

if (dn_p_max>=dn_m_max):
	dn=dn_p_max
else:
	dn=dn_m_max

print("Dn is ",dn)

ks_stat=(n**0.5+0.12+0.11/(n)**0.5)*dn_m_max

print("The K-S statistics is",ks_stat)

### ending K-S test

### starting serial test

d=3

lcg_rv1_serial=[]
for i in range(n*d):
    z_i=(a*z_0+c)%m
    lcg_rv1_serial.append(z_i/m)
    #print(z_i/m)
    z_0=z_i

serial_u = [[0 for i in range(d)] for j in range(n)]

for i in range(n):
	serial_u[i][0]=mcg_rv1_serial[d*i]
	serial_u[i][1]=mcg_rv1_serial[d*i+1]
	serial_u[i][2]=mcg_rv1_serial[d*i+2]

#print(serial_u)

# use 2d_list to store Nj values

# for j in range(k):
# 	for num in range(n_rep):
# 		if (rv_storer[num]>=(j-1)/k) and (rv_storer[num]<(j/k)):
# 			k_list[j]=k_list[j]+1

# # create a k list and initialize all k intervals to 0

k=5

N_j_d=[[[0 for i in range(k)] for j in range(k)] for m in range(k)]

print(len(N_j_d))

for j0 in range(k):
	for j1 in range(k):
		for j2 in range(k):
			for i in range(n):
				if serial_u[i][0]>=(j0-1)/k and serial_u[i][0]<(j0/k) and serial_u[i][1]>=(j1-1)/k and serial_u[i][1]<(j1/k) and serial_u[i][2]>=(j2-1)/k and serial_u[i][2]<(j2/k):
					N_j_d[j0][j1][j2]+=1

#print(N_j_d)

# calculate chi-stat

chi_d=0

for i in range(k):
	for j in range(k):
		for m in range(k):
			chi_d+=(N_j_d[i][j][m]-n/(k**d))

chi_d=k**d/n

print("the serial test statistics is ", chi_d)

## ending serial test

## starting runs test

a_matrix=[[4529.4,9044.9,13568,18091,22615,27892],[9044.9,18097,27139,36187,45234,55789],[13568,27139,40721,54281,67852,83685],[18091,36187,54281,72414,90470,111580],[22615,45234,67852,90470,113262,139476],[27892,55789,83685,111580,139476,172860]]

b_mat=[1/6,5/24,11/120,19/720,29/540,1/840]

n=5000

r=[0,0,0,0,0,n]

Rn=0


lcg_rv1_runs=[]
for i in range(n):
    z_i=(a*z_0+c)%m
    lcg_rv1_runs.append(z_i/m)
    #print(z_i/m)
    z_0=z_i

### mrg

# n=5000

for i in range(3,n):
    z1_i=(a1_1*z1_0+a1_2*z1_1+a1_3*z1_2)%(2**32-209)
    z2_i=(a2_1*z2_0+a2_2*z2_1+a2_3*z2_2)%(2**32-22853)
    z1_list.append(z1_i)
    z2_list.append(z2_i)

mcg_rv1_auto=[]
for i in range(n):
    y_i=(z1_list[i]-z2_list[i])%(2**32-209)
    if y_i>0:
        u_i=y_i/(2**32-208)
    else:
        u_i=(2**32-209)/(2**32-208)
    mcg_rv1_serial.append(u_i)



for i in range(6):
	for j in range(6):
		sum=a_matrix[i][j]*(r[i]-n*b_mat[i])*(r[j]-n*b_mat[j])
		Rn+=sum

Rn=1/n*Rn

print("Rn is ",Rn)

## starting autocorrelation test

lag=6
n=1000

h=math.floor(((n-1)/lag)-1)

p_j=0

# for i in range(h):
# 	sum_auto=lcg_rv1[i*lag]*lcg_rv1[(i+1)*lag]
# 	p_j+=sum_auto

for i in range(h):
	sum_auto=mcg_rv1[i*lag]*mcg_rv1[(i+1)*lag]
	p_j+=sum_auto

p_j=12/(1+h)*p_j-3

var_pj=(13*h+7)/((h+1)**2)

Aj=p_j/math.sqrt(var_pj)

print("Aj is", Aj)


### Buffon's needle estimation

B_n=0

L=0.5

n=1000


# boardWidth = 50
# needleLength = 50*L
# numberOfNeedles = n
    
# myPen = turtle.Turtle()
# myPen.hideturtle()
# myPen.speed(0)

# y=180
# #Draw floor boards
# for i in range(0,10):
#   myPen.penup()
#   myPen.goto(-200,y)
#   myPen.pendown()
#   myPen.goto(200,y)
#   y-=boardWidth

# #Draw Needles
# myPen.color("#f442d1")
# for needle in range(0,numberOfNeedles):
#   x=mcg_rv1[needle]*360-180
#   y=mcg_rv2[needle]*360-180
#   angle=mcg_rv3[needle]*360
#   myPen.penup()
#   myPen.goto(x,y)
#   myPen.setheading(angle)
#   myPen.pendown()
#   myPen.forward(needleLength)
  
# print("L = " + str(needleLength))
# print("N = " + str(numberOfNeedles))
# print("W = " + str(boardWidth))

### Dart throwing estimation

# pi, via dart throwing

totalThrows=1000     # how many times we throw the dart
throwsInsideCircle = 0 # starting value for a counter
for throw in range(totalThrows): # loop
  x = mcg_rv1[throw]*2 -1 # a random x value between [-1,1]
  y = mcg_rv2[throw]*2 -1 # random y value
  if(x*x + y*y <= 1.0): # if inside the circle (distance squared)
    throwsInsideCircle += 1 

pi = (4.0*throwsInsideCircle)/totalThrows
print("pi is",pi) 
