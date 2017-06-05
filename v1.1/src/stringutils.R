# Wrappers functions for common string operations. 
# 
# Author: doreper
###############################################################################

stringutils = new.env(hash=T)

stringutils$getFilesWithPattern = function(directory, pattern)
{
    return(file.path(directory, list.files(directory, pattern))) 
}

stringutils$getFilesWithExtension = function(directory, extension)
{
#	if(substr(extension,1,1)!=".")
#	{
#		stop("expects extension starting with .")
#	}
	pattern = paste0("\\", extension, "$")
	return(file.path(directory, list.files(directory, pattern))) 
}

stringutils$getbaseprefix=function(fname, postfix)
{
	split = strsplit(basename(fname), postfix)
	prefix = unlist(lapply(split, "[", 1))
}
