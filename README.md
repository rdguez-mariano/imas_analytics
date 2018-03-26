# IMAS analytics

A compilation of MATLAB functions that were used to apply, analyse and visualise Image Matching by Affine Simulation (**IMAS**) methods based on available state-of-the-art Scale Invariant Image Matching (**SIIM**) methods.

These functions were conceived to analyse results from the following articles:
- [Covering the Space of Tilts](https://rdguez-mariano.github.io/pages/imas) - SIIMS
- [Fast affine invariant image matching](https://rdguez-mariano.github.io/pages/hyperdescriptors) - IPOL
- [Affine invariant image comparison under repetitive structures](https://rdguez-mariano.github.io/pages/acdesc) - ICIP

## Getting Started

The first requirement is that a working IMAS C++ code should be available in the folder "IMAS_cpp" or, optionally, a working executable in the folder "siim-matcher". All functions were intended to be used with the first option (C++ code) but you may do small modifications to get them working with the second one. Normally, the first time that you run any function that requires compilation through MEX will be compiled automatically and hopefully you won't need go into details.

### Prerequisites

If you want to use the [USAC](http://www.cs.unc.edu/~rraguram/usac/) algorithm you'll need to install some libraries first. For example, on ubuntu you should do something like:
```bash
sudo apt-get update
sudo apt-get install libblas-dev liblapack-dev libconfig-dev libconfig++-dev
```
The fact of using USAC is controlled by a compilation variable inside **"perform_IMAS.m"** called `USAC`.

Also, if you want to use local descriptors provided by OPENCV then you need to download and compile opencv v3.2.0 and opencv_contrib v3.2.0. Then modify the function **"perform_IMAS.m"** so as paths to folder libraries and modules do follow the proper path.

## Computing near optimal coverings as in the [IPOL article](https://rdguez-mariano.github.io/pages/hyperdescriptors)
Inside the folder **"opt_covering"** you can find the C++ code for finding near optimal coverings based on radius (initial visibility) and the region to cover. This code doesn't have any documentation available yet but there exist a pseudo-code that has been published in IPOL [(See its web page)](https://rdguez-mariano.github.io/pages/hyperdescriptors).
### Getting started
Compilation in linux :
```bash
cd opt_covering && mkdir -p build && cd build && cmake .. && make
```
If you have problems to compile with the png library please set the `PNG` variable to OFF in the **"CMakeLists.txt"**.
### Running example
Executing the program is very simple, the following example in bash is clear enough:
 ```bash
 r=1.6           # radius of disks ( 1 <= r < Inf )
 region=4.5      # region to be covered ( r < region < Inf )

 N=2             # number of groups of concentric disks
                 # (1 <= N << Inf)

 epsilon=1.025   # Parameter for tuning annulus discretisation.
                 # Used to check for covering conditions
                 # ( 1 < epsilon < Inf )

 ./opt_covering $r $region $epsilon $N
 ```

 In the case of a radius equal to 1.4 you might want to add an extra group of concentric disks for the same region.
 ```bash
 ./opt_covering 1.4 4.5 1.025 3
 ```

## Running, Visualising and Analysing examples

### Force to re-compile perform_IMAS
If you have modified some files in the IMAS_cpp folder, you may want to compile and execute in a single call setting the compile option to true.
```matlab
opts_IMAS.compile = true;
data_matches = perform_IMAS(im1,im2,opts_IMAS);
```

### Sanity check of found formulas for disks boundaries in the [SIIMS article](https://rdguez-mariano.github.io/pages/imas)
Let A be a fixed element of the space of tilts `$ A = T_4 R_{\frac{pi}{4}} $`. Every random C (`$ C = T_2 R_\theta $`) applied to A (i.e. CA) must lie in the boundary of the disk of center A and radius 2. The decomposition of CA is achieved numerically through the SVD decomposition and independent of our formulas.

```matlab
optionsdraw1.drawtitle = false;
optionsdraw1.tilt = 4*sqrt(2);
theta = pi/4;
draw_ball_pol((4),{theta},log((2)),optionsdraw1);
R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
A=[4 0;0 1]*R;

theta = 2*pi/rand; % pi/7; et  pi/8;
T=[2 0;0 1];
R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
C = T*R;


opts_comp.human_output = true;
opts_comp.epsilon = 10^(-10);

[ t_sim,theta_sim,~,~ ] = affine_decomp(C*A , opts_comp);
plot(log(t_sim).*cos(theta_sim),log(t_sim).*sin(theta_sim),'.b','MarkerSize',12);
```

### Numerical composition of elements in the space of tilts
The following code shows what it seems to be a cyclic Group in the Space of tilts. The generator of the group in this case is `$T_2 R_{\frac{pi}{7}}$`.

```matlab
t = 2;
theta = pi/7;

opts_comp.human_output = true;
opts_comp.epsilon = 10^(-10);
comp = 30;

T=[2 0;0 1];
R = [cos(theta) sin(theta); -sin(theta) cos(theta)];

tvec = [t];
psicell = {[theta]};

temp = T*R;
for n=1:comp % new levels to add
    temp = T*R*temp;
    [ t_sim,theta_sim,~,~ ] = affine_decomp(temp , opts_comp);
    theta = [tvec t_sim];
    psicell = [psicell theta_sim];
    tvec = [tvec t_sim];
end

optionsdraw2.tilt = [4*sqrt(2)];
optionsdraw2.interactive = false;
optionsdraw2.draw_circle_eachtilt = false;
optionsdraw2.draw_regions = true;
optionsdraw2.draw_regions_tilt = 4*sqrt(2);

figure;
[tvec1,psicell1] = draw_ball_pol(tvec,psicell,log(sqrt(2)),optionsdraw2);

figure;
draw_ball3d(tvec,psicell,log(sqrt(2)),optionsdraw2);
```


### FAIR-SURF example

First, load two gray images into `im1` and `im2`. For example, do
```matlab
im1 = double(imread('./image_BD/adam1.png'));
im2 = double(imread('./image_BD/adam2.png'));
```

Then get the exact covering proposed for FAIR-SURF by typing,
```matlab
[ tvec, psicell, radius, region ] = get_literature_covering('FAIR-SURF fixed tilts covering');
```

or, in general, if you want to test a near optimal covering you could type instead,
```matlab
radius = 1.6;
[ tvec, psicell, region ] = get_feasible_covering( radius );
```

Compute the area ratio and prepare options to be passed.
```matlab
val = 0; count =0;
for i=1:length(tvec)
    t=tvec(i);
    numphi=length(psicell{i});
    count = count + numphi;
    val = val + numphi/t;
end

opts_IMAS.tvec1 = tvec;
opts_IMAS.tvec2 = tvec;
opts_IMAS.psicell1 = psicell;
opts_IMAS.psicell2 = psicell;
opts_IMAS.desc_type = 2;
opts_IMAS.match_ratio = 0.8;

optionsdraw.tilt = [region];
optionsdraw.draw_regions = true;
optionsdraw.draw_regions_tilt = max(tvec(:))*radius;
```

Visualise the corresponding covering in 2D and 3D.
```matlab
figure;
draw_ball_pol(tvec,psicell,log(((radius))),optionsdraw);
htitle = get(gca,'Title');
title( [get(htitle,'String') ' / Area ratio = ' num2str(val) ' / Tilt \leq ' num2str(region)]);

figure;
draw_ball3d(tvec,psicell,log(radius),optionsdraw);
title( ['simulations = ' num2str(count) ' / radius = ' num2str(radius) ' / Area ratio = ' num2str(val) ' / Tilt \leq ' num2str(region)]);
```

Execute FAIR-SURF on our two images.
```matlab
data_matches = perform_IMAS(im1,im2,opts_IMAS);
```

Visualise in a zenith view from what pairs of simulated images the matches were coming.
```matlab
opts_zenith.threshold = 15;
opts_zenith.linewidth = 2;
opts_zenith.markersize = 15;
plot_IMAS_zenith( data_matches, opts_zenith );
```


### Tilt tolerance test applied in the [SIIMS article](https://rdguez-mariano.github.io/pages/imas)
This test is based on the assumption that ORSA Homography never fails in two ways: first, if it exists an underlying homography that explains the matches then ORSA will find it; and if ORSA tells that an homography explains the matches then that homography is truly the underlying one. Also, another assumption is that planes related by the underlying homography are not composed of repetitive structures.

[See this article on HAL for more details](https://hal.archives-ouvertes.fr/hal-01589522)

Execute this test by typing :

```matlab
opts_maxtilt.draw_on = true;
anglevec = 0:10:170;
tvec = 1.2:0.1:2.4;

figure;
[Mbool, Mval, RADIUS ] =   tolerance_tests_SIIMS(tvec,anglevec,opts_maxtilt);

title(['Transition Tilt Tolerance t_{max} = ' num2str(RADIUS) ' )']);
```

### Tilt tolerance test with repetitive structures applied in the [ICIP article](https://rdguez-mariano.github.io/pages/acdesc)

This test is based on the assumption that repetitive structures are present on both images. The goal is to keep track on the number and ratio of true matches.

[See this article](/)

Execute this test by typing :

```matlab
num_test = 100;


im1 = double(imread('./image_BD/enami.png'));
im2 = double(imread('./image_BD/porta.png'));
im3 = double(imread('./image_BD/adam1.png'));

matrixgood = [];
matrixratio = [];
tvec = 1:0.1:3;
for t=tvec
    ratiovec = 0;
    goodvec = 0;
    badvec = 0;
    grafratio = [];
    grafgood = [];
    grafbad = [];
    figure;
    for i=1:num_test
        [goodm, badm] = tolerance_test_ICIP( im1, im2, im3, 100, 3, t, randi([0 179]) );
        ratio = goodm./(goodm+badm);
        ratio(isnan(ratio)) = 0.5; % definition 0/0 = 0.5

        ratiovec = ratiovec + ratio;
        goodvec = goodvec + goodm;
        badvec = badvec +badm;

        grafratio = [grafratio; ratio];
        grafgood = [grafgood; goodm];
        grafbad = [grafbad; badm];

        subplot(2,2,1);
        plot(1:i,grafratio(:,1),'-xc');hold on;
        plot(1:i,grafratio(:,6),':xc');
        plot(1:i,grafratio(:,2),'-+m');
        plot(1:i,grafratio(:,7),':+m');
        plot(1:i,grafratio(:,3),'-sb');
        plot(1:i,grafratio(:,4),'-ok');
        %plot(1:i,grafratio(:,5),'-*r');hold off;
        legend('SIFT L1 0.8','SIFT L1 0.6','RootSIFT 0.8', 'RootSIFT 0.6', 'Weights', 'Quantised 0.3');%,'SURF');
        title('Ratio')
        axis auto;

        subplot(2,2,2);
        plot(1:i,cumsum(grafratio(:,1),1)./((1:i)'),'-xc');hold on;
        plot(1:i,cumsum(grafratio(:,6),1)./((1:i)'),':xc');
        plot(1:i,cumsum(grafratio(:,2),1)./((1:i)'),'-+m');
        plot(1:i,cumsum(grafratio(:,7),1)./((1:i)'),':+m');
        plot(1:i,cumsum(grafratio(:,3),1)./((1:i)'),'-sb');
        plot(1:i,cumsum(grafratio(:,4),1)./((1:i)'),'-ok');
        % plot(1:i,cumsum(grafratio(:,5),1)./((1:i)'),'-*r');hold off;
        legend('SIFT L1 0.8','SIFT L1 0.6','RootSIFT 0.8', 'RootSIFT 0.6', 'Weights', 'Quantised 0.3');%,'SURF');
        title('Mean Ratio')
        axis auto;

        subplot(2,2,3);
        plot(1:i,grafgood(:,1),'-xc');hold on;
        plot(1:i,grafgood(:,6),':xc');
        plot(1:i,grafgood(:,2),'-+m');
        plot(1:i,grafgood(:,7),':+m');
        plot(1:i,grafgood(:,3),'-sb');
        plot(1:i,grafgood(:,4),'-ok');
        %plot(1:i,grafgood(:,5),'-*r');hold off;
        legend('SIFT L1 0.8','SIFT L1 0.6','RootSIFT 0.8', 'RootSIFT 0.6', 'Weights', 'Quantised 0.3');%,'SURF');
        title('Good')
        axis auto;

        subplot(2,2,4);
        plot(1:i,cumsum(grafgood(:,1),1)./((1:i)'),'-xc');hold on;
        plot(1:i,cumsum(grafgood(:,6),1)./((1:i)'),':xc');
        plot(1:i,cumsum(grafgood(:,2),1)./((1:i)'),'-+m');
        plot(1:i,cumsum(grafgood(:,7),1)./((1:i)'),':+m');
        plot(1:i,cumsum(grafgood(:,3),1)./((1:i)'),'-sb');
        plot(1:i,cumsum(grafgood(:,4),1)./((1:i)'),'-ok');
        % plot(1:i,cumsum(grafgood(:,5),1)./((1:i)'),'-*r');hold off;
        legend('SIFT L1 0.8','SIFT L1 0.6','RootSIFT 0.8', 'RootSIFT 0.6', 'Weights', 'Quantised 0.3');%,'SURF');
        title('Mean good')
        axis auto;


        refresh;
        pause(0.2);
    end

    %means
    matrixratio = [matrixratio; ratiovec/num_test];
    matrixgood = [matrixgood; goodvec/num_test];
    badvec/num_test

    if (false)
        save(['matlab_t_' num2str(t) '.mat']);
    end
end


figure;
plot(tvec,matrixratio(:,1),'-xc');hold on;
plot(tvec,matrixratio(:,6),':xc');
plot(tvec,matrixratio(:,2),'-+m');
plot(tvec,matrixratio(:,7),':+m');
plot(tvec,matrixratio(:,3),'-sb');
plot(tvec,matrixratio(:,4),'-ok');
        legend('SIFT L1 0.8','SIFT L1 0.6','RootSIFT 0.8', 'RootSIFT 0.6', 'AC-W', 'AC-Q 0.3');%,'SURF');
        xlabel('Viewpoint angle')
        ylabel('Mean');
        xticks(1:0.1:3)

        label = {};
        for tilt=1:0.1:3
            label = [label {[num2str(round(180*acos(1/tilt)/pi)) '\circ']}];    
        end
        xticklabels(label)
title('Ratio')

figure;
plot(tvec,matrixgood(:,1),'-xc');hold on;
plot(tvec,matrixgood(:,6),':xc');
plot(tvec,matrixgood(:,2),'-+m');
plot(tvec,matrixgood(:,7),':+m');
plot(tvec,matrixgood(:,3),'-sb');
plot(tvec,matrixgood(:,4),'-ok');
        legend('SIFT L1 0.8','SIFT L1 0.6','RootSIFT 0.8', 'RootSIFT 0.6', 'AC-W', 'AC-Q 0.3');%,'SURF');
        xlabel('Viewpoint angle')
        ylabel('Mean');
        xticks(1:0.1:3)

        label = {};
        for tilt=1:0.1:3
            label = [label {[num2str(round(180*acos(1/tilt)/pi)) '\circ']}];    
        end
        xticklabels(label)
title('Good')
```

### Similarity test

```matlab
im1 = double(imread('./image_BD/adam1.png'));
im2 = double(imread('./image_BD/enami.png'));


similarity_test( im1, im2, 100, 3, 2, 45 )
```

## Author

[Mariano Rodríguez](https://rdguez-mariano.github.io/)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Some of these functions were also conceived with:
  * Julie Delon
  * Jean-Michel Morel
  * Rafael Grompone von Gioi
* The folder simu_tilts_cpp were taken and modified accordingly from the ASIFT project developed by [Guoshen Yu](http://www.cmap.polytechnique.fr/~yu/).
