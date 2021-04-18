using MAT
using ImageFiltering
using ImageView, Images
using ImageContrastAdjustment
bandsfile = matread("matbands.mat")

function create_test(image, s)
    img1 = imfilter(image, Kernel.gaussian(1/s))
    img2 = imfilter(img1, (1/(s^2)).*ones(s, s))
    imout = img2[1:s:end, 1:s:end]
    return imout
end

bands = bandsfile["bands"]

B2 = bands["B2"]
B3 = bands["B3"]
B4 = bands["B4"]
B4 = bands["B4"]
B5 = bands["B5"]
B6 = bands["B6"]
B7 = bands["B7"]
B8 = bands["B8"]
B8a = bands["B8a"]
B9 = bands["B9"]
B11 = bands["B11"]
B12 = bands["B12"]

#B2_down = bandsfile["B2_down"]
#B3_down = bandsfile["B3_down"]
#B4_down = bandsfile["B4_down"]
#B8_down = bandsfile["B8_down"]

#RGB_og = cat(B4, B3, B2, dims=3)
#imshow(RGB_og)

B2_down = create_test(B2, 2)
B3_down = create_test(B3, 2)
B4_down = create_test(B4, 2)
B5_down = create_test(B5, 2)
B6_down = create_test(B6, 2)
B7_down = create_test(B7, 2)
B8_down = create_test(B8, 2)
B8a_down = create_test(B8a, 2)
B9_down = create_test(B9, 2)
B11_down = create_test(B11, 2)
B12_down = create_test(B12, 2)

#=
b2 = adjust_histogram(B2_down, LinearStretching(dst_minval = 0, dst_maxval = 1))
b3 = adjust_histogram(B3_down, LinearStretching(dst_minval = 0, dst_maxval = 1))
b4 = adjust_histogram(B4_down, LinearStretching(dst_minval = 0, dst_maxval = 1))

RGB_down = cat(b2, b3, b4, dims=3)
colors = colorview(RGB, RGB.(RGB_down))
imshow(RGB_down)
=#

#write(bands_out, "B2_down", B2_down, "B3_down", B3_down, "B4_down", B4_down, "B8_down", B8_down)

matwrite("downsampled.mat", Dict("B2_down" => B2_down, "B3_down" => B3_down, "B4_down" => B4_down, "B5_down" => B5_down, "B6_down" => B6_down, "B7_down" => B7_down, "B8_down" => B8_down, "B8a_down" => B8a_down, "B9_down" => B9_down, "B11_down" => B11_down, "B12_down" => B12_down))
