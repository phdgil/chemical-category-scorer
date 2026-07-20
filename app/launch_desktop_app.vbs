Dim shell, fso, appDir, launcher, localPy314, localPy313, localPy312, localPy311
Set shell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
appDir = fso.GetParentFolderName(WScript.ScriptFullName)
launcher = Chr(34) & appDir & "\launch_desktop_app.pyw" & Chr(34)
shell.CurrentDirectory = appDir

localPy314 = shell.ExpandEnvironmentStrings("%LocalAppData%") & "\Programs\Python\Python314\pythonw.exe"
localPy313 = shell.ExpandEnvironmentStrings("%LocalAppData%") & "\Programs\Python\Python313\pythonw.exe"
localPy312 = shell.ExpandEnvironmentStrings("%LocalAppData%") & "\Programs\Python\Python312\pythonw.exe"
localPy311 = shell.ExpandEnvironmentStrings("%LocalAppData%") & "\Programs\Python\Python311\pythonw.exe"

If fso.FileExists(localPy314) Then
    shell.Run Chr(34) & localPy314 & Chr(34) & " " & launcher, 0, False
    WScript.Quit 0
End If
If fso.FileExists(localPy313) Then
    shell.Run Chr(34) & localPy313 & Chr(34) & " " & launcher, 0, False
    WScript.Quit 0
End If
If fso.FileExists(localPy312) Then
    shell.Run Chr(34) & localPy312 & Chr(34) & " " & launcher, 0, False
    WScript.Quit 0
End If
If fso.FileExists(localPy311) Then
    shell.Run Chr(34) & localPy311 & Chr(34) & " " & launcher, 0, False
    WScript.Quit 0
End If

If fso.FileExists("C:\Python314\pythonw.exe") Then
    shell.Run Chr(34) & "C:\Python314\pythonw.exe" & Chr(34) & " " & launcher, 0, False
    WScript.Quit 0
End If
If fso.FileExists("C:\Python313\pythonw.exe") Then
    shell.Run Chr(34) & "C:\Python313\pythonw.exe" & Chr(34) & " " & launcher, 0, False
    WScript.Quit 0
End If
If fso.FileExists("C:\Python312\pythonw.exe") Then
    shell.Run Chr(34) & "C:\Python312\pythonw.exe" & Chr(34) & " " & launcher, 0, False
    WScript.Quit 0
End If
If fso.FileExists("C:\Python311\pythonw.exe") Then
    shell.Run Chr(34) & "C:\Python311\pythonw.exe" & Chr(34) & " " & launcher, 0, False
    WScript.Quit 0
End If

shell.Run "pyw -3 " & launcher, 0, False
