from tsvd import  *
import  numpy as np
import  time
import  matplotlib.pyplot as plt
import h5py
from scipy.io import  loadmat
from scipy.io import  savemat
#loading spectrum data, size I * J * K * T

rho = 0.15
input_filename = '../Xtrue_rho'+str(rho)+'.mat'
data = loadmat(input_filename)

Spectrum_tens = np.array(data['Xacc'])
Spectrum_tens = Spectrum_tens.transpose([0,2,1]) # Transpose to I*T*J*K
#for 3D data, re-arrange the data to size M*T*K




n1, n2, n3 = Spectrum_tens.shape #load data size

#slab by slab analysis
Xcube = Spectrum_tens

X = Tensor(Xcube)
Xfrob = tfrobnorm(X)

tic = time.time()
U,S,V = tsvd(X,full=False)
toc = time.time()

S = np.diag(S.array()[:,:,0])
tElapsed = toc - tic
print('Time elapsed for tSVD: {:.3f}'.format(tElapsed))
plt.figure(figsize=(10,5))
plt.semilogy(S)
# plt.xticks(np.arange(0, min(n1,n2), step=20))
plt.rcParams.update({'font.size': 15})
plt.show()

total_sum = np.sum(S)
pwrs = np.cumsum(S) / total_sum
# idx = np.where(pwrs < 0.99)
idx = np.where(pwrs<0.99)
# print('75% power rank: '+ np.str(idx[0][-1] + 1))
print('99% power rank:'+str(idx[0][-1]+1))

tubalrank = int(idx[0][-1]+1)

np.random.seed(0)

# tube = False
# rho = 0.8 #data missing rate, i.e., sampling ratio = 1 - rho
#
# if (tube is False):
#     mask = np.random.rand(n1, n2, n3)
#     mask[mask > rho] = 1
#     mask[mask <= rho] = 0
#     mask = mask.astype(int)
# else:
#     mask = np.random.rand(n1, n2)
#     mask[mask > rho] = 1
#     mask[mask <= rho] = 0
#     mask = mask.astype(int)
#     mask = np.repeat(mask[:, :, np.newaxis], n3, axis=2)

# mask mat
data = loadmat(('../mask_rho'+str(rho)+'.mat'))
mask = np.array(data['Wacc'])
mask = mask.transpose([0,2,1])



#Tensor SVD
################################ TOUCAN #######################################
rank = 1
# fun = lambda L: [0, tfrobnorm(L - X) / Xfrob]
fun = lambda Q,k: [0, tfrobnorm_array(Q.array()[:,k,:] - X.array()[:,k,:]) / tfrobnorm_array(X.array()[:,k,:])]
# fun = lambda Q,k : [0,0]
Y_hat_assem = np.zeros([n1,n2,n3])


Xcube_freq = Spectrum_tens[:, :, :]








# inputmat
data = loadmat('../Y_rho'+str(rho)+'.mat')
Y = np.array(data['Yacc'])#Y = Tensor(Xcube_freq * mask)
Y = Y.transpose([0,2,1])
Y = Tensor(Y) # M*T*K
Y_hat_toucan, U_toucan,stats_toucan, tElapsed_toucan = toucan(Y,mask,rank,tube=False,outer=1,cgtol=1e-6,mode='online',fun=fun,
                                                 randomOrder=False,verbose=False)
Y_hat_assem = Y_hat_toucan.array()

toucan_err_curve = np.zeros(n2)

for time_slot in range(n2):
    Y_hat_cube_t = Y_hat_assem[:, time_slot,:]
    Ground_truth_t = Spectrum_tens[:, time_slot,:]
    toucan_err_curve[time_slot] = tfrobnorm_array(Y_hat_cube_t - Ground_truth_t) / tfrobnorm_array(Ground_truth_t)

print('TOUCAN Time: {:4f}'.format(tElapsed_toucan))


mdic = {"NMSEcomp":toucan_err_curve}
savemat(("../../result/NMSE_TOUCAN_rho"+str(rho)+".mat"),mdic)


times_toucan = stats_toucan[1:,-1]

toucan_err = tfrobnorm(Y_hat_toucan - X) / Xfrob
print('NRMSE TOUCAN: {:6f}'.format(toucan_err))
plt.semilogy(toucan_err_curve)
plt.title('TOUCAN: Recovered Tensor NRMSE')
plt.xlabel('Iteration')
plt.show()
# plt.scatter(np.arange(0,len(cgiter_toucan[1:])),cgiter_toucan[1:])
# plt.show()

