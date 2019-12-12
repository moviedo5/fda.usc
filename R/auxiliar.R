  ind.list=function (lista, ind)
{
    len <- length(lista)
    for (i in 1:len) lista[[i]] <- lista[[i]][ind, , drop = FALSE]
    lista
}
