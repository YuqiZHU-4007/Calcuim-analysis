function callback_readme()
open('F:\DUlab\FC analyse\Ca_img_processingº”…Ÿ÷°∞Ê\readme.txt')
end
function opentxt(filename)
   [~, name, ext] = fileparts(filename); 
   fprintf('You have requested file: %s\n', [name ext]);

   if exist(filename, 'file') == 2
     fprintf('Opening in MATLAB Editor: %s\n', [name ext]);
     edit(filename);
   else
      wh = which(filename);
      if ~isempty(wh)
         fprintf('Opening in MATLAB Editor: %s\n', wh);
         edit(wh);
      else
        warning('MATLAB:fileNotFound', ...
                'File was not found: %s', [name ext]);
      end
   end
end