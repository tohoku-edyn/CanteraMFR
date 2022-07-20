# CanteraMFR

Simulate MFR using Cantera.

## Ubuntu 18.04 or 20.04

### Prepare software libraries
Install fmt

    sudo apt install libfmt-dev

Install boost

    sudo apt install libboost-dev
    
Install git

    sudo apt install git    

Install Cantera

    sudo apt install software-properties-common
    sudo apt-add-repository ppa:cantera-team/cantera
    sudo apt install cantera-python3 cantera-dev
    
### Compile CanteraMFR

    git clone https://github.com/Youhichka/CanteraMFR.git
    cd CanteraMFR
    mkdir build
    cd build
    cmake ../
    make

**Note: If your Ubuntu version is 18.04, don't forget to modify 'CMakeLists.txt'.**
    
### Run an example of CanteraMFR for Ammonia
    cd ..
    mkdir work
    cp -r examples/ammonia/ work/
    cd work
    ./mfr

### Set conditions
Please modify "inputs.yaml".

# FAQ on the MFR system
- [Computational model for weak flame in MFR](http://www.ifs.tohoku.ac.jp/enerdyn/en/research/mfr-faq1.html)
- [Temperature difference between weak flame and wall](http://www.ifs.tohoku.ac.jp/enerdyn/en/research/mfr-faq2.html)
- [Negligibly small impact of radical quenching](http://www.ifs.tohoku.ac.jp/enerdyn/en/research/mfr-faq3.html)

# [How to cite](../../wiki/How-to-cite)


