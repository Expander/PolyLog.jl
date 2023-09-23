import DelimitedFiles

# read file row-wise into pairs of complex numbers
function read_from(filename, T=Float64)
    data = DelimitedFiles.readdlm(filename, T)
    z1 = map(z -> Complex(z[1], z[2]), zip(data[:,1], data[:,2]))
    z2 = map(z -> Complex(z[1], z[2]), zip(data[:,3], data[:,4]))
    hcat(z1, z2)
end
