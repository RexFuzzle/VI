function saveplot(filename,fignum,H,W)
if exist('filename')==0
filename=input('Please type filename without extension: ','s');
end
if exist('fignum')==0
fignum=input('Which figure number: ');
end
if exist('H')==0
H=input('Height: ');
end
if exist('W')==0
W=input('Width: ');
end
set(figure(fignum),'PaperUnits','inches')
set(figure(fignum),'PaperOrientation','portrait');
set(figure(fignum),'PaperSize',[H,W])
set(figure(fignum),'PaperPosition',[0,0,W,H])
print(figure(fignum),'-dtex',sprintf('%s.tex',filename));
