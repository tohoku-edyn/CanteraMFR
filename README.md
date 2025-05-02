# CanteraMFR

Simulate MFR using Cantera.

## Ubuntu 22.04

### Prepare software libraries
Update & Upgrade

    sudo apt update
    sudo apt upgrade
    
Install fmt

    sudo apt install libfmt-dev

Install boost

    sudo apt install libboost-dev
    
Install git & gh

    sudo apt install git 
    sudo apt install gh
    gh auth login

Install Cantera

    sudo apt install software-properties-common
    sudo apt-add-repository ppa:cantera-team/cantera
    sudo apt install cantera-python3 cantera-dev
    
### Compile CanteraMFR

    gh repo clone tohoku-edyn/CanteraMFR
    cd CanteraMFR
    mkdir build
    cd build
    cmake ../
    make

**Note: If your Ubuntu version is 18.04 or 20.04, don't forget to modify 'CMakeLists.txt'.**
    
### Run an example of CanteraMFR for Ammonia
    cd ..
    cp -r examples/ammonia2024/* work/
    cd work
    ./mfr inputs.yaml

### Set conditions
Please modify "inputs.yaml".

### For Windows users
The method shown below is one example. There might be other methods.

    1. Install "Windows Subsystem for Linux" from Microsoft Store.
    2. Install "Ubuntu" from Microsoft Store.
    3. Enable "Virtual Machine Platform" in Windows Control Panel.
    4. Run Ubuntu.
    5. Install Cantera.
         sudo apt update
         sudo apt install python3-pip
         python3 -m pip install cantera
    6. Follow instructions for Ubuntu above.

# FAQ on the MFR system
- [Computational model for weak flame in MFR](http://www.ifs.tohoku.ac.jp/enerdyn/en/research/mfr-faq1.html)
- [Temperature difference between weak flame and wall](http://www.ifs.tohoku.ac.jp/enerdyn/en/research/mfr-faq2.html)
- [Negligibly small impact of radical quenching](http://www.ifs.tohoku.ac.jp/enerdyn/en/research/mfr-faq3.html)

# [How to cite](../../wiki/How-to-cite)


