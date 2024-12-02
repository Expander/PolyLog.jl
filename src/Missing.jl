for f in (:li0, :li1, :li2, :li3, :li4, :li5, :li6, :reli1, :reli2, :reli3, :reli4)
    @eval $(f)(::Missing) = missing
end

for f in (:li, :reli)
    @eval $(f)(::Integer, ::Missing) = missing
end
