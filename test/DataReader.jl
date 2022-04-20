import DelimitedFiles

function read_from(filename)
    data = DelimitedFiles.readdlm(filename)
    input  = map(z -> Complex(z[1], z[2]), zip(data[:,1], data[:,2]))
    output = map(z -> Complex(z[1], z[2]), zip(data[:,3], data[:,4]))
    hcat(input, output)
end
