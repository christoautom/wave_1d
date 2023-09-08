# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 14:38:25 2023

@author: Christophe Roman
"""


import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
#mpl.use("pgf")

N=199
index=np.arange(0,N)
Dt=.001
t_spam=np.arange(0,15,Dt)


u_=list()
du_=list()
a_=list()
dx_=list()
for i in index :
    exec("u_{}=sp.symbols('u_{}')".format(i,i))
    eval("u_.append(u_{})".format(i))
    exec("du_{}=sp.symbols('du_{}')".format(i,i))
    eval("du_.append(du_{})".format(i))
    exec("a_{}=sp.symbols('a_{}')".format(i,i))
    exec("a_.append(a_{})".format(i,i))
    exec("dx_{}=sp.symbols('dx_{}')".format(i,i))
    exec("dx_.append(dx_{})".format(i,i))

g_1=sp.symbols('g_1')
b_1=sp.symbols('b_1')
p_=[g_1,b_1]
L1_=list()
L2_=list()
L3_=list()
L4_=list()
L5_=0
for i in np.arange(1,N-1) :
    eval('L1_.append(.5*du_{}**2*dx_{})'.format(i,i))
    eval('L2_.append(-a_{}/2/6*((u_{}-u_{})/dx_{})**2*dx_{})'.format(i-1,i,i-1,i,i))
    eval('L3_.append(-a_{}/2*4/6*((u_{}-u_{})/(dx_{}+dx_{}))**2*dx_{})'.format(i,i+1,i-1,i,i+1,i))
    eval('L4_.append(-a_{}/2/6*((u_{}-u_{})/dx_{})**2*dx_{})'.format(i+1,i+1,i,i+1,i))
#exec("L5_=[.5*a_{}*du_{}+.5*a_{}*du_{}] ".format(N-1,N-1,0,0))
exec("L5_=[.5*a_{}*du_{}**2/g_1+.5*a_{}*du_{}**2/b_1] ".format(N-1,N-1,0,0))

#exec("L5_=[.5*du_{}**2/g_1+.5*du_{}**2/b_1] ".format(N-1,0))
L=np.sum(np.array(L1_)+np.array(L2_)+np.array(L3_)+np.array(L4_))+np.array(L5_)

E=list()
dE=list()
A_=np.zeros([N,N],dtype='object')
EA=np.zeros([1,N],dtype='object')
for i in index :
    eval("E.append(sp.diff(L[0],u_{}))".format(i))
    eval("dE.append(sp.diff(L[0],du_{}))".format(i))
    for j in index :
        exec("A_[{},{}]=(sp.diff(E[{}],u_{}))".format(i,j,i,j))
    exec("EA[0,{}]=(sp.diff(dE[{}],du_{}))".format(i,i,i))
f_A=sp.Lambda((tuple(a_+dx_+p_)),(1/EA*A_.T).T)
f_E=sp.Lambda((tuple(a_+dx_+p_)),EA)
#print_latex(Matrix(A))
def f_R(q):
    R=np.diag(q)
    return R

def f_B(N):
    B=np.zeros([N,1])
    B[N-1,0]=1
    return B




def spam_system(a,q,p,f,Dx,t_spam,ut_ref,K,u0,u1) :
    Dt_spam=np.diff(t_spam)
    u=np.empty([N,len(t_spam)])
    dot_u=np.empty([N,len(t_spam)])
    U=np.empty([1,len(t_spam)])
    eta=np.empty([1,len(t_spam)])
    A=(np.array(f_A(*(a+Dx+p)))).astype('float')
    E=(np.array(f_E(*(a+Dx+p)))).astype('float')
    R=f_R(q)/E
    B=f_B(N)/E.T
    F=np.array(f)/E
    k=0
    u[:,0]=u0
    dot_u[:,0]=u1
    I=np.eye(N)
    eta[0,0]=0
    U[0,k]=-K[1]*eta[0,k]-K[0]* (dot_u[N-1,k]-ut_ref)
    for t in Dt_spam:
        k=k+1
        dot_u[:,k]=np.linalg.inv(I+t*R)@(dot_u[:,k-1]+t*A@u[:,k-1]+(t*B.T*U[0,k-1])[0]+t*F[0])
        u[:,k]=u[:,k-1]+t*dot_u[:,k]       
        eta[0,k]=eta[0,k-1]+t*(dot_u[N-1,k-1]-ut_ref)
        U[0,k]=-K[1]*eta[0,k]-K[0]* (dot_u[N-1,k]-ut_ref)
    return dot_u,u,U,eta




    
if 1:
    u0=np.zeros([N])
    u1=np.zeros([N])
    Dx=list(index*0+1/N)
    xx=np.cumsum(Dx)
    a=list(np.sin(xx*2)+2)
    f=list(np.sin(xx*2*3.14))
    f[0]=1
    f[N-1]=-1
    q=list(.01+.1*xx*2)
    q[0]=1
    q[N-1]=1
    p=[20,20]
    ut_ref=5
    K=[20,10]


if 0:
    u0=np.zeros([N])
    u1=np.zeros([N])
    Dx=list(index*0+1/N)
    xx=np.cumsum(Dx)
    a=list(np.sin(xx*2)+2)
    f=list(np.sin(xx*2*3.14))
    f[0]=1
    f[N-1]=-1
    q=list(.01+.1*xx*2)
    q[0]=1
    q[N-1]=1
    p=[20,20]
    ut_ref=5
    K=[0,0]

dot_u,u,U,eta=spam_system(a,q,p,f,Dx,t_spam,ut_ref,K,u0,u1)

print('SIMULATION DONE')    
x_spam=np.arange(0,N)/(N-1)
ts,x=np.meshgrid(t_spam, x_spam)
if 0:

    
   

    vertices = np.column_stack((ts.flatten(), x.flatten(), dot_u.flatten() ))
    faces = []
    num_rows, num_cols = x.shape
    for i in range(num_rows - 1):
        for j in range(num_cols - 1):
            v1 = i * num_cols + j
            v2 = v1 + 1
            v3 = v1 + num_cols
            v4 = v2 + num_cols
            faces.append([v1, v2, v3])
            faces.append([v2, v4, v3])
    
 
    faces = np.array(faces)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_trisurf(vertices[:, 0], vertices[:, 1], vertices[:, 2], triangles=faces, cmap='viridis')

    ax.set_xlabel(r'time $t$')
    ax.set_ylabel(r'space $x$')
    ax.set_zlabel(r'$v_t(x,t)$')  
    ax.set_box_aspect([2, 1, 1])
    fig.tight_layout(pad=.5)
    #eval('''plt.savefig('{}',dpi=300, bbox_inches='tight')'''.format(ff))
    plt.savefig('fig_ut.png',dpi=300, bbox_inches='tight') 
    print('print u_t done')
if 0:   
    vertices = np.column_stack((ts.flatten(), x.flatten(), u.flatten() ))
    faces = []
    num_rows, num_cols = x.shape
    for i in range(num_rows - 1):
        for j in range(num_cols - 1):
            v1 = i * num_cols + j
            v2 = v1 + 1
            v3 = v1 + num_cols
            v4 = v2 + num_cols
            faces.append([v1, v2, v3])
            faces.append([v2, v4, v3])
    
 
    faces = np.array(faces)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_trisurf(vertices[:, 0], vertices[:, 1], vertices[:, 2], triangles=faces, cmap='viridis')

    ax.set_xlabel(r'time $t$')
    ax.set_ylabel(r'space $x$')
    ax.set_zlabel(r'$v(x,t)$')  
    ax.set_box_aspect([2, 1, 1])
    fig.tight_layout(pad=.5)
    plt.savefig('fig_u.png',dpi=300, bbox_inches='tight') 
    print('print u done')
if 0:   
    vertices = np.column_stack((ts.flatten(), x.flatten(), np.gradient(u,axis=0).flatten() ))
    faces = []
    num_rows, num_cols = x.shape
    for i in range(num_rows - 1):
        for j in range(num_cols - 1):
            v1 = i * num_cols + j
            v2 = v1 + 1
            v3 = v1 + num_cols
            v4 = v2 + num_cols
            faces.append([v1, v2, v3])
            faces.append([v2, v4, v3])
    
 
    faces = np.array(faces)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_trisurf(vertices[:, 0], vertices[:, 1], vertices[:, 2], triangles=faces, cmap='viridis')

    ax.set_xlabel(r'time $t$')
    ax.set_ylabel(r'space $x$')
    ax.set_zlabel(r'$v_x(x,t)$')  
    ax.set_box_aspect([2, 1, 1])
    fig.tight_layout(pad=.5)
    plt.savefig('fig_ux.png',dpi=300, bbox_inches='tight') 
    print('print u_x done')
if 1:
    
    fig, ax = plt.subplots(figsize=(4.5, 2.5))
    line1,=ax.plot(t_spam,U[0,:],linewidth=.5,label=r'U(t)')
     

    ax.set_xlabel(r'time $t$')
    ax.set_ylabel(r'U(t)')
    ax.legend(handles=[line1])
    fig.tight_layout(pad=.5)  
    plt.savefig('fig1_U.pdf', backend='pgf')
    
    fig, ax = plt.subplots(figsize=(4.5, 2.5))
    line1,=ax.plot(t_spam,dot_u[0,:],linewidth=.5,linestyle='dashed',label=r'$v_t(0,t)$')
    line2,=ax.plot(t_spam,dot_u[N-1,:],linewidth=.5,label=r'$v_t(1,t)$')
    ax.set_xlabel(r'time $t$')
    ax.set_ylabel(r'velocity')
    ax.legend(handles=[line1, line2])
    fig.tight_layout(pad=.5)   
    plt.savefig('fig1_B.pdf', backend='pgf')
    
    obj=dot_u-ut_ref
    fig, ax = plt.subplots(figsize=(4.5, 2.5))
    line1,=ax.plot(t_spam,obj[0,:],linewidth=1,linestyle='dashed',label=r'$v_t(0,t)-v_t^{ref}$')
    line2,=ax.plot(t_spam,obj[N-1,:],linewidth=1,label=r'$v_t(1,t)-v_t^{ref}$')
    xobj=np.zeros_like(t_spam)
    i=-1
    for t in t_spam:
        i=i+1
        xobj[i]=np.sum(Dx*obj[:,i])
    line3,=ax.plot(t_spam,xobj,linewidth=1,linestyle='dotted',label=r'$\int_0^x (v_t(x,t)-v_t^{ref})dx$')
    ax.set_xlabel(r'time $t$')
    ax.set_ylabel(r'objective')
    ax.legend(handles=[line1, line2, line3])
    fig.tight_layout(pad=.5)   
    plt.savefig('fig1_obj.pdf', backend='pgf')
    # plt.figure()
    # plt.plot(f)
    # plt.plot(q)
    # plt.plot(a)


    # # Plot the surface
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # ax.plot_surface(t, x, dot_u, vmin=dot_u.min() * 2, cmap=cm.Blues)
    # ax.set(xticklabels=[],
    #    yticklabels=[],
    #    zticklabels=[])

   
    # fig.tight_layout(pad=.5)
    
    # plt.savefig('fig1.pdf', backend='pgf')
