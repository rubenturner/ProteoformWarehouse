#!/bin/bash


# Some initial dependencies
xcode-select --install
pip3 install requests --user
pip3 install pyqt5

# Downloading homebrew installer
curl -fsSL -o install.sh https://raw.githubusercontent.com/Homebrew/install/master/install.sh
/bin/bash install.sh

# Configuring pip installer
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py

# Install more dependencies
brew install pyqt5
sudo pip3 install requests --upgrade
pip3 install PySide2
