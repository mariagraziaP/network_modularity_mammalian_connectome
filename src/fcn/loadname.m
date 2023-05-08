function newvar=loadname(fullfilename,name)

var=load(fullfilename);
structname=who('-file',fullfilename);
aux=sprintf('var.%s',structname{1});
newvar=eval(aux);