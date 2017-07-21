@everywhere using PTstartup
@everywhere println(getMyGamma())

# function fWriteBinary(stream::IOStream, param::Float64, M::Int64)
#     for i = 1:M
#         write(stream, i + param)
#         # flush(stream)
#         write(stream, 3.14)
#     end
#     return
# end

# function fWriteBinary32(stream::IOStream, param::Float32, M::Int32)
#     for i = 1:M
#         write(stream, [i + param, 3.14])
#         # flush(stream)
#         # write(stream, 3.14)
#     end
#     return
# end


# function fWriteText(stream::IOStream, param::Float64, M::Int64)
#     for i = 1:M
#         println(stream, i + param)
#         println(stream, 3.14)
#         flush(stream)
#     end
#     return
# end

# testfile1 = "testfile1.dat"
# testfile2 = "testfile2.dat"
# stream1 = open(testfile1, "w")
# stream2 = open(testfile2, "w")
# M = 9
# M32 = Int32(15000000)
# param = 0.5
# param32 = Float32(0.5)

# @time fWriteBinary(stream1, param, M)
# # @time fWriteBinary32(stream2, param32, M32)
# # @time fWriteText(stream2, param, M)

# function fRead(fname)
#     outfile_content = open(fname, "r")
#     i = 0

#     while !eof(outfile_content)
#         i += 1
#         println(read(outfile_content, Float64, (1, 2)))
#     end
#     return i
# end


# @time counter = fRead("/home/wikash/Documents/msc-thesis_project/Code/20170525-write-files/testfile1.dat")
# println(counter)


#  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# a = zeros(Float64, 2, M);
# for i = 1:4
# for i = 1:3
#     # println(read(outfile_content, Float64, 1)[1])
#     # println(a[i])
    # a[:, i] = read(outfile_content, Float64, (1, 2))
# end
# print(2 * q)
# print(2 * w)
# print(r)
# print(a[:, 1])
# print(a[:, 2])
# print(a[:, 4])
# print(a[:, 1123213])

# g = read(outfile_content, Int64, 1)
# h = read(outfile_content, Int64, (1,2))
# i = read(outfile_content, Int64, 1)
# j = read(outfile_content, Int64, 1)
# println(g)
# println(h)
# println(i)
# println(j)



# outfile = "outfile.dat"
# # writing to files is very similar:
# f = open(outfile, "w")
# # both print and println can be used as usual but with f as their first arugment
# a = 3
# b = [4, 3]
# c = [4.0, 3.0, 2]
# d = 44
# write(f, a)
# write(f, b)
# write(f, d)
# # println(f, a)
# # print(f, "more content")
# # print(f, " more on the same line")
# # flush(f)
# # println(f, b)
# # flush(f)
# # write(f, c)
# close(f)

# we can then check the content of the file written
# "do" above just creates an anonymous function and passes it to open
# we can use the same logic to pass readall and thereby succinctly
# open, read and close a file in one line
# outfile_content = open(readall, outfile, "r")
# println(repr(outfile_content))
#> "some content\nmore content more on the same line"
# outfile_content = open(outfile, "r")
# g = read(outfile_content, Int64, 1)
# h = read(outfile_content, Int64, (1,2))
# i = read(outfile_content, Int64, 1)
# # j = read(outfile_content, Int64, 1)
# println(g)
# println(h)
# println(i)
# println(j)
# close(outfile)
