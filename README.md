# HEMat_old
Implementation of homormophic matrix multiplication based on (Huang &amp; Zong, 2023) and (Zhu et al. 2023) in HElib 

[2023a] Huang, H., Zong, H.: Secure matrix multiplication based on fully homomorphic encryption. The Journal of Supercomputing 79(5), 5064–5085 (2023).
[2023b] Zhu, L., Hua, Q., Chen, Y., Jin, H.: Secure outsourced matrix multiplication with fully homomorphic encryption. In: European Symposium on Research in Computer Security, pp. 249–269. Springer, Cham (2023). 

## Step 1. Building and installing HElib
### General prerequisites
* Ubuntu 22.04 LTS
* GNU make >= 4.3
* g++ >= 11.3.0
* cmake >= 3.22
* pthreads
* git >=2.36
### Instructions
1. Installing basis tools:
```
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install patchelf
sudo apt-get install m4
sudo apt-get install git
sudo apt-get install cmake
```
2. Clone HElib from gitbub:
```
sudo git clone https://github.com/homenc/HElib.git
```
3. Install HElib
```
cd HElib
sudo mkdir build
cd build
sudo cmake -DPACKAGE_BUILD=ON -DCMAKE_INSTALL_PREFIX=/home/usr/helib_install ..
sudo make -j16 -lpthread
sudo make install
```
Here HElib is installed in `/home/usr/helib_install`. 



## Step 2. Download the code from `HEMat_old` and build

### Instructions

1. Clone the code
```
git clone https://github.com/XiaopengZheng/HEMat_old.git
```

2. Build the code and run (HElib has been installed in `/home/usr/helib_install`)
```
cd HEMat
sudo mkdir build
cd build
sudo cmake -DPACKAGE_BUILD=ON -DCMAKE_INSTALL_PREFIX=/home/usr/helib_install ..
sudo make
cd bin
./main
```




