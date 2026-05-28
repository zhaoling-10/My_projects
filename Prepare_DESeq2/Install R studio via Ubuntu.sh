############

Install R

############

Step 1) Install R in Ubuntu

sudo apt update
sudo apt install r-base r-base-dev -y

Step 2) Check version

R --version

Step 3) Install RStudio Desktop in Ubuntu

wget https://download1.rstudio.org/electron/jammy/amd64/rstudio-2023.12.1-402-amd64.deb

sudo apt install ./rstudio-2023.12.1-402-amd64.deb


Step 4) SUID sandbox permission

sudo chown root:root /usr/lib/rstudio/chrome-sandbox
sudo chmod 4755 /usr/lib/rstudio/chrome-sandbox


Step 5) Launch RStudio

rstudio

Step 6) OPTIONAL (not mandatory): Fix GPU rendering issues in RStudio

6.1) disable GPU for RStudio

rstudio --disable-gpu

6.2) Make it permanent (desktop launcher)

mkdir -p ~/.local/share/applications
cp /usr/share/applications/rstudio.desktop ~/.local/share/applications/
sed -i 's/^Exec=rstudio.*/Exec=rstudio --disable-gpu %F/' ~/.local/share/applications/rstudio.desktop
update-desktop-database ~/.local/share/applications
