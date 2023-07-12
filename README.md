# CpHMD_v3
Git commit of the new CpHMD version
This repository contains the old versions of CpHMD_v2 (normal, verlet, US, US_verlet) and now has the new aditions for CpHMD_v3. 

Also contains the Sts and Delphitools which will be completed in the future! 

### -- LOGBOOK -- ###
## 12/07/2023 ##
- Introduction of AMBER FF (by J. Sequeira) into CpHMD v3 framework *(Debbug still on going, use carefully)*
    - Creation of Delphi_tools for AMBER
- Integration of PybindE tool (by J. Vitorino) into the CpHMD Framework.*(Debbug still on going, use carefully)*


## 05/07/2023 ##
- Introduction of dendrimer portion of the code
- Added dendrimer blocks to delphi, sts and FF GROMOS54a7
- corrected C-termini when capped to 8 due to the adition of SPB

## 04/07/2023 ##
- Added block for cobimetinib and cobimetinib-EP to delphi, sts and FF GROMOS54a7
- Virtual site code was added to CpHMD_v3, allowing it do deal with EPs and others.
- Corrected a bug with the sed that corrected O1 on C-terminal
