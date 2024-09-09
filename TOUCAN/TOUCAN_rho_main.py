from tsvd import  *
import  numpy as np
import  time
import  matplotlib.pyplot as plt
import h5py
from scipy.io import  loadmat
from scipy.io import  savemat
#loading spectrum data, size I * J * K * T
file_path = 'G:\PapersAndSlides\MyPaper2024\MyPaper2024\TSP_2024\TSP_Recode\spectrum_data'
file_name = '\Param_R8_sigma8_K64_T600_v0.1_var_TOUCAN.mat'
input_filename = file_path + file_name
data = loadmat(input_filename)

Spectrum_tens = np.array(data['X4DT'])
Spectrum_tens = Spectrum_tens.transpose([0,3,1,2]) # Transpose to I*T*J*K



n1, n2, n3, n4 = Spectrum_tens.shape

#slab by slab analysis
Xcube = Spectrum_tens[:, :, :, 15]

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
print('99% power rank:'+np.str(idx[0][-1]+1))

tubalrank = int(idx[0][-1]+1)

np.random.seed(0)

tube = False
rho = np.linspace(start=0.9,stop=0.8,num=11) #data missing rate, i.e., sampling ratio = 1 - rho
length = np.size(rho)






#Tensor SVD
################################ TOUCAN #######################################
rank = 1
toucan_err_curve = np.zeros((length,n2))
for pp in range(11):
    p = rho[pp]
    if (tube is False):
        mask = np.random.rand(n1, n2, n3)
        mask[mask > p] = 1
        mask[mask <= p] = 0
        mask = mask.astype(int)
    else:
        mask = np.random.rand(n1, n2)
        mask[mask > p] = 1
        mask[mask <= p] = 0
        mask = mask.astype(int)
        mask = np.repeat(mask[:, :, np.newaxis], n3, axis=2)

    sig = 0

# fun = lambda L: [0, tfrobnorm(L - X) / Xfrob]
    fun = lambda Q,k: [0, tfrobnorm_array(Q.array()[:,k,:] - X.array()[:,k,:]) / tfrobnorm_array(X.array()[:,k,:])]
    # fun = lambda Q,k : [0,0]
    Y_hat_assem = np.zeros([n1,n2,n3,n4])

    for freq in range(n4):
        Xcube_freq = Spectrum_tens[:, :, :, freq]
        Y = Tensor(Xcube_freq * mask)
        Y_hat_toucan, U_toucan,stats_toucan, tElapsed_toucan = toucan(Y,mask,rank,tube=False,outer=1,cgtol=1e-6,mode='online',fun=fun,
                                                         randomOrder=False,verbose=False)
        Y_hat_assem[:,:,:,freq] = Y_hat_toucan.array()

    for time_slot in range(n2):
        Y_hat_cube_t = Y_hat_assem[:, time_slot,:,:]
        Ground_truth_t = Spectrum_tens[:, time_slot,:,:]
        toucan_err_curve[pp,time_slot] = tfrobnorm_array(Y_hat_cube_t - Ground_truth_t) / tfrobnorm_array(Ground_truth_t)

    print('TOUCAN Time: {:4f}'.format(tElapsed_toucan))
    mdic = {"NMSE_TOUCAN": toucan_err_curve}
    savemat("G:\PapersAndSlides\MyPaper2024\MyPaper2024\TSP_2024\TSP_Recode\Results\VaryingR\R8T600rho0.1-0.2v0.1_var.mat",
            mdic)


# times_toucan = stats_toucan[1:,-1]
#
# toucan_err = tfrobnorm(Y_hat_toucan - X) / Xfrob
# print('NRMSE TOUCAN: {:6f}'.format(toucan_err))
# plt.semilogy(toucan_err_curve)
# plt.title('TOUCAN: Recovered Tensor NRMSE')
# plt.xlabel('Iteration')
# plt.show()
# plt.scatter(np.arange(0,len(cgiter_toucan[1:])),cgiter_toucan[1:])
# plt.show()

