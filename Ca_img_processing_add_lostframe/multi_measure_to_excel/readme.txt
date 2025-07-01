跑 ROI_multimeasure_to_excel_ZYQ.ijm
1、打开imageJ；导入文件；run
2、选择配准后的文件夹（每层子文件夹的上一级）
3、选择excel表格输出路径
4、选择ROI set所在文件夹
结果：在输出路径输出同名excel.xls
注意：
每层roi set的名字与配准后的文件夹同名；
roi set所在文件夹中‘*.zip’都是ROI set；
配准后的文件夹只包含配准后的tiff文件（建议检查文件数是否一致）

跑Excel_xls2xlsx.cls
1、在excel输出文件夹创建新的excel表；
2、打开开发工具；visual basic；导入文件Excel_xls2xlsx
3、run
结果：在问价路径下创建文件夹xlsx；包含各层的表格.xlsx

跑getRoi_from_xls_in_diff_file.m
1、选择xlsx所在文件夹
2、选择输出mat所在文件夹
结果：在输出路径输出‘activities_new.mat’
注意：
1、activities_new.mat的顺序和输入路径excel表一致，如果每层不是按顺序的，先排好
2、选择哪一行，如果multimeasure的输出变了要调整

%20181129,ZYQ