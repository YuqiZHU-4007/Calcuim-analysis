VERSION 1.0 CLASS
BEGIN
  MultiUse = -1  'True
END
Attribute VB_Name = "ThisWorkbook"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = True
Sub xls2xlsx()
Dim FilePath, MyFile, iPath, Name, OutPath As String
iPath = ThisWorkbook.Path
OutPath = Dir(iPath & "\xlsx", vbDirectory)
If OutPath = "" Then
    MkDir (iPath & "\xlsx")
End If
MyFile = Dir(iPath & "\*.xls")

If MyFile <> "" Then
Do
    On Error Resume Next
    If MyFile = ThisWorkbook.Name Then MyFile = Dir
    Workbooks.Open (iPath & "\" & MyFile)
    MyFile = Replace(MyFile, ".xls", ".xlsx")
    Name = "\" & MyFile
    FilePath = iPath & "\xlsx" & Name
    Application.ScreenUpdating = False
    ActiveWorkbook.SaveAs Filename:=FilePath, FileFormat:=xlOpenXMLWorkbook, CreateBackup:=False
    Workbooks(MyFile).Close True
    Application.ScreenUpdating = True
    MyFile = Dir
Loop While MyFile <> ""
End If
End Sub
