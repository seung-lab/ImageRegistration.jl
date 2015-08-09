# Thomas Macrina
# 150804
# Working from https://github.com/dfdx/PiecewiseAffineTransforms.jl/blob/master/examples/ex.jl

println("Simple example of transforming person's face from one shape to another")
Pkg.update()
Pkg.add("VoronoiDelaunay")
Pkg.add("PiecewiseAffineTransforms")

using Images
using ImageView
using MAT
using VoronoiDelaunay
using PiecewiseAffineTransforms
using Color
using FixedPointNumbers

const IMG_HEIGHT = 480
pwd()

img_path = joinpath(pwd(), "image.jpg")

# some auxilary functions for reading data
rawdata(img) = convert(Array{Float64, 3}, data(separate(img)))
xy2ij(shape, height) = [height .- shape[:, 2] shape[:, 1]]

read_image(path) = rawdata(imread(path))
read_shape(path) = xy2ij(matread(path)["annotations"], IMG_HEIGHT)

# read 2 images: src_img will be transformed to resemble dst_img
src_img = read_image(img_path)
# dst_img = convert(Array, separate(imread(joinpath(data_dir, dst_img_name))))
dst_img = read_image(img_path)

a =   [291.0995  347.8602
  301.7107  366.5858
  325.4298  372.8277
  355.3908  375.3244
  373.4922  362.8407
  421.5546  357.8472
  435.2867  369.7068
  464.6235  370.9551
  488.9668  362.8407
  505.8199  341.6183
  312.3218  327.8862
  322.3088  337.8732
  337.9135  341.6183
  353.5182  337.2490
  367.2503  326.0137
  354.7666  324.7653
  337.9135  324.1411
  322.9330  326.0137
  430.2932  324.7653
  441.5286  331.6313
  459.0058  336.0006
  472.1138  331.0072
  485.8459  319.1476
  474.6105  317.2750
  461.5026  316.6508
  446.5221  319.1476
  419.0579  331.6313
  419.0579  294.8043
  435.9109  281.0722
  440.9044  242.3726
  424.0514  254.2321
  404.7016  238.6274
  380.9824  256.7289
  365.3778  239.2516
  363.5052  276.0787
  382.8550  296.0527
  384.1034  332.8797
  331.6717  199.3036
  384.7276  204.9213
  398.4597  204.2971
  416.5611  204.2971
  459.6300  194.9343
  421.5546  164.3492
  404.7016  159.3557
  378.4857  161.8524
  515.1827  197.4311
  510.1892  174.3361
  498.9538  153.1138
  487.7185  132.5156
  464.6235  105.0514
  443.4012   90.6951
  402.8290   85.0774
  352.2698   91.3192
  327.3023  109.4207
  308.5767  133.1398
  292.3479  156.8589
  282.9850  183.6990
  277.3674  211.7874];
src_shape = xy2ij(a, IMG_HEIGHT);

b = [  278.3889  354.2037
  296.9074  366.7963
  313.9444  369.7593
  350.2407  369.7593
  366.5370  354.2037
  419.1296  345.3148
  434.6852  360.1296
  465.0556  355.6852
  484.3148  351.2407
  503.5741  329.7593
  299.1296  326.0556
  308.7593  334.2037
  324.3148  338.6482
  340.6111  334.2037
  356.1667  321.6111
  345.0556  317.9074
  325.7963  316.4259
  310.9815  320.1296
  428.0185  314.9444
  439.1296  320.8704
  458.3889  324.5741
  472.4630  317.9074
  481.3519  305.3148
  469.5000  300.1296
  454.6852  302.3519
  441.3519  305.3148
  412.4630  314.9444
  413.9444  274.9444
  425.7963  254.2037
  425.7963  219.3889
  410.9815  231.2407
  393.2037  214.9444
  368.0185  237.9074
  355.4259  219.3889
  353.9444  254.2037
  373.9444  274.2037
  376.1667  316.4259
  322.8333  160.1296
  366.5370  177.9074
  385.7963  171.9815
  404.3148  174.9444
  437.6481  152.7222
  395.4259  140.1296
  379.1296  141.6111
  364.3148  143.0926
  500.6111  173.4630
  489.5000  146.0556
  477.6481  120.1296
  459.8704   91.9815
  440.6111   69.0185
  413.2037   59.3889
  383.5741   57.9074
  345.7963   60.8704
  314.6852   73.4630
  291.7222   96.4259
  275.4259  124.5741
  259.8704  156.4259
  250.2407  203.0926];
dst_shape = xy2ij(b, IMG_HEIGHT);

# get triangulation
trigs = delaunayindexes(src_shape)

println("Original images and shapes (press Enter to continue)")
triplot(src_img, src_shape, trigs)
triplot(dst_img, dst_shape, trigs)
readline(STDIN)

println("Simple warp (press Enter to continue)")
@time warped_simple = pa_warp(src_img, src_shape, dst_shape, trigs)
view(src_img)
view(warped_simple)
readline(STDIN)

println("Preparing warp parameters")
@time pa_params = pa_warp_params(dst_shape, trigs, (480, 640))

println("Prepared warp (press Enter to finish)")
@time warped_prepared = pa_warp(pa_params, src_img, src_shape)
view(src_img)
view(warped_prepared)
readline(STDIN)

function gettriangles