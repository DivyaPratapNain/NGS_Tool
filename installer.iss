[Setup]
AppName=Diff. Exp. Tool
AppVersion=1.0
DefaultDirName={pf}\DiffExpTool
DefaultGroupName=DiffExpTool
OutputBaseFilename=DiffExpToolInstaller
Compression=lzma2
SolidCompression=yes
PrivilegesRequired=admin

[Files]
; Include your main app executable
Source: "dist\app.exe"; DestDir: "{app}"; Flags: ignoreversion
; Include the dependency installer script
Source: "install_deps.py"; DestDir: "{app}"; Flags: ignoreversion
; Include the batch file for running WSL
Source: "run_wsl.bat"; DestDir: "{app}"; Flags: ignoreversion

[Run]
; Run the batch file for WSL execution
Filename: "{app}\run_wsl.bat"; Parameters: ""; Flags: runmaximized waituntilterminated

[Icons]
; Create a desktop shortcut
Name: "{commondesktop}\Diff. Exp. Tool"; Filename: "{app}\app.exe"
; Create a start menu shortcut
Name: "{group}\Diff. Exp. Tool"; Filename: "{app}\app.exe"
