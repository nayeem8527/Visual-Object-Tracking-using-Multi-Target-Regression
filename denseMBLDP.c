
/*

denseMBLDP extract Histograms of Multi-Block LDP on a regular dense grid and computed on image I after color projection


Usage
------

[dmbldp , infodmbldp] =  densembldp(I , [options] );


Inputs
-------

I                              Input image (ny x nx x [3]) in UINT8 format. 
options 

       scale                   Multi-Scale vector (1 x nscale) (default scale = 1) where scale(i) = s is the size's factor to apply to each 9 blocks
                               in the LDP computation, i = 1,...,nscale

	   deltax                  Division step in the x-axis for the grid (default deltax = floor(nx*min(scale))) 
	   deltay                  Division step in the y-axis for the grid (default deltay = floor(ny*min(scale)))

       n                       Derivative order of ldp (n > 1, default n = 2)

       ldporient               Orientation for ldp computation (1 x norient), where norien <= 4 and ldporient(i) = {0,1,2,3}
	                           ldporient(i) = 0 <=> (I(x_1+1,y_1)-I(x_1,y_1))*(I(x_0+1,y_0)-I(x_0,y_0))< 0
	                           ldporient(i) = 1 <=> (I(x_1+1,y_1+1)-I(x_1,y_1))*(I(x_0+1,y_0+1)-I(x_0,y_0))< 0
	                           ldporient(i) = 2 <=> (I(x_1,y_1+1)-I(x_1,y_1))*(I(x_0,y_0+1)-I(x_0,y_0))< 0
	                           ldporient(i) = 3 <=> (I(x_1-1,y_1+1)-I(x_1,y_1))*(I(x_0-1,y_0+1)-I(x_0,y_0))< 0
	                           (default ldporient = [0 , 1 , 2 , 3])

       color                   0 : force gray-scale (dimcolor = 1, default), 1 : RGB-LDP (dimcolor = 3), 2 : nRGB-LDP (dimcolor = 3), 3 : Opponent-LDP (dimcolor = 3), 
                               4 : nOpponent-LDP (dimcolor = 2), 5 : Hue-LDP (dimcolor = 1)

	   patchsize               Size of the patch where the descriptor is computed (default patchsize = 16 )	  


       improvedLDP             0 for default 8 bits encoding, 1 for 9 bits encoding (8 + central area)

       kernelx                 Kernel in x-direction for computing the gradient (default kernelx = [1])
       kernely                 Kernel in y-direction for computing the gradient (default kernely = [1])
                               If numel(kernelx) ~1 || numel(kernely) ~1, compute the magnitude, i.e. sqrt(dx*dx + dy*dy)
                               For example if kernelx = [1 2 1 ; 0 0 0 ; -1 -2 -1] and kernely = kernelx', this is equivalent to compute filter
                               the image with a sobel filter.


       rmextremebins           Force to zero bin = {0 , {255,58,9}} if  rmextremebins = 1 (default rmextremebins = 1)

       norm                    Normalization : norm = 0 <=> no normalization, norm = 1 <=> v=v/(sum(v)+epsi), norm = 2 <=> v=v/sqrt(sum(v²)+epsi²), 
                               norm = 3 <=> v=sqrt(v/(sum(v)+epsi)) , norm = 3 <=> L2-clamped (default norm = 1)

       clamp                   Clamping value (default clamp = 0.2)


Outputs
-------

dmbldp                         MBLDP patches descriptors (nbin x nb_pts) nb_pts = deltax*deltay*nscale*dimcolor 
                               in double format where nbin = ([256,59,10]*(improvedLDP+1)) (depending on improvedLDP) 
							   

infodmbldp                     MBLDP patches informations(4 x nb_pts) where nb_pts = deltax*deltay*nscale*dimcolor
                               infodmbldp(1,i) = y
						       infodmbldp(2,i) = x
							   infodmbldp(3,i) = scale
							   infodmbldp(4,i) = color    
							   infodmbldp(5,i) = nyscale                 
							   infodmbldp(6,i) = nxscale                


To compile
----------


mex  -g -output denseMBLDP.dll denseMBLDP.c

mex  -f mexopts_intel10.bat -output denseMBLDP.dll denseMBLDP.c




Example 1
---------

I                     = imread('02769_Right_StudentOffice.jpeg');
%I                     = imread('image_0174.jpg');

options.scale         = [1 , 2];
options.deltax        = 10;
options.deltay        = 10;
options.n             = 3;
options.ldporient     = [0 , 1 , 2 , 3];
options.color         = 0;
options.patchsize     = 16;
options.improvedLDP   = 0;
options.rmextremebins = 1;
options.norm          = 4;
options.clamp         = 0.2;


[dmbldp , infodmbldp] =  denseMBLDP(I , [options] );

figure(1)
imagesc(dmbldp)


figure(2)
imagesc(I)
colormap(gray)
hold on
plot(infodmbldp(2 , :) , infodmbldp(1 , :) , 'r+')
hold off


Example 2
---------

I                     = imread('image_0174.jpg');

options.scale         = [1,2,3];
options.deltax        = 12;
options.deltay        = 12;
options.n             = 2;
options.ldporient     = [0 , 1 , 2 , 3];
options.color         = 0;
options.patchsize     = 20;
options.improvedLDP   = 0;
options.rmextremebins = 1;
options.norm          = 4;
options.clamp         = 0.2;


[dmbldp , infodmbldp] =  denseMBLDP(I , options );

half                  = options.patchsize/2;

xr                    = [infodmbldp(2, :)-half ; infodmbldp(2, :)-half ; infodmbldp(2, :)+ half ; infodmbldp(2, :)+ half ; infodmbldp(2, :)-half] + 1.5;
yr                    = [infodmbldp(1, :)-half ; infodmbldp(1, :)+half ; infodmbldp(1, :)+ half ; infodmbldp(1, :)- half ; infodmbldp(1, :)-half] + 1.5;



figure(1)
imagesc(I)
colormap(gray)
hold on
plot(infodmbldp(2 , :)+1.5 , infodmbldp(1 , :)+1.5 , 'r+')
plot(xr , yr , 'b')
hold off


figure(2)
imagesc(dmbldp)




Author : Sébastien PARIS : sebastien.paris@lsis.org
-------  Date : 09/06/2010

Reference : [1] Baochang Zhang, Yongsheng Gao, Sanqiang Zhao, and Jianzhuang Liu, 
---------       Local Derivative Pattern Versus Local BinaryPattern: Face Recognition With High-Order Local Pattern Descriptor
                IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 19, NO. 2, FEBRUARY 2010
           

*/


#include <time.h>
#include <math.h>
#include <mex.h>


#define tiny 1e-8
#define verytiny 1e-15
#define sqrt3 1.73205080756888
#define invsqrt2 0.707106781186547
#define invsqrt3 0.577350269189626
#define invsqrt6 0.408248290463863

#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif


struct opts
{
	double *scale;
	int nscale;
	double *kernelx;
	int    kyx;
	int    kxx;
	double *kernely;
	int    kyy;
	int    kxy;
	int    deltax;
	int    deltay;
	int    n;
    int    *ldporient;
	int    norient;
	int    color;
	int    patchsize;
	int    improvedLDP;
	int    rmextremebins;
	int    norm;
	double clamp;
};


/*-------------------------------------------------------------------------------------------------------------- */
/* Function prototypes */

void rgb2gray(unsigned char * , int , int , double *);
void rgb2nrgb(unsigned char * , int , int , double *);
void rgb2opponent(unsigned char * , int , int , double *);
void rgb2nopponent(unsigned char * , int , int , double *);
void rgb2hue(unsigned char * , int , int , double *);


void MakeIntegralImage(double *, double *, int , int , double *);
double ld(double * , int , int , int , int , int , int );
double Area(double * , int , int , int , int , int );
void conv2(double *, double *, int  , int , int  , int , int , int , double * , double * );
void compute_hmbldp(double * , int , int , int , int , int , int ,  int , double *, int , struct opts);
void denseMBLDP(double * , int  , int , struct opts , double * , double *  , int , int);

/*-------------------------------------------------------------------------------------------------------------- */

#ifdef MATLAB_MEX_FILE
void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{  
	unsigned char *im;
	const int *dimim;
	double *I;
	double *feat , *des;

	double	kernelx_default[1]   = {1};
	double	kernely_default[1]   = {1};	
	int     ldporient_default[4] = {0,1,2,3};

	int i , ny , nx , nynx , dimcolor = 1  , npts ,  npt_ker = 1 , npt_ldp = 4, Nbins = 256;
	struct opts options;
	mxArray *mxtemp;
	double *tmp , temp , scalemin = 1.0;
	int tempint , dummy;

	options.nscale        = 1;
	options.kyx           = 1;
	options.kxx           = 1;
	options.kxy           = 1;
	options.kyy           = 1;
	options.deltax        = 10;
	options.deltay        = 10;
	options.n             = 2;
	options.norient       = 4;
	options.color         = 0;
	options.patchsize     = 16;
	options.improvedLDP   = 0;
	options.rmextremebins = 1;
	options.norm          = 0;
	options.clamp         = 0.2;


	/* Input 1  */


	if( (mxGetNumberOfDimensions(prhs[0]) > 3) || !mxIsUint8(prhs[0]) )
	{
		mexErrMsgTxt("I must be (ny x nx x [3]) in UINT8 format\n");
	}


	im          = (unsigned char *)mxGetData(prhs[0]);
	dimim       = mxGetDimensions(prhs[0]);
	ny          = dimim[0];
	nx          = dimim[1];
	nynx        = ny*nx;

	/* Input 2  */


	if ((nrhs > 1) && !mxIsEmpty(prhs[1]) )
	{

		mxtemp                            = mxGetField(prhs[1] , 0 , "scale");
		if(mxtemp != NULL)
		{
			if( mxGetM(mxtemp) != 1 )
			{
				mexErrMsgTxt("spyr must be (1 x nscale) in double format\n");
			}
			options.scale                 = mxGetPr(mxtemp);
			options.nscale                = mxGetN(mxtemp);
		}

		else
		{
			options.nscale                = 1;
			options.scale                 = (double *)mxMalloc(sizeof(double));
			options.scale[0]              = 1.0;
		}

		mxtemp                            = mxGetField( prhs[1], 0, "kernelx" );
		if(mxtemp != NULL)
		{

			options.kernelx               = mxGetPr(mxtemp);
			options.kyx                   = (int) mxGetM(mxtemp);
			options.kxx                   = (int) mxGetN(mxtemp);

		}

		else
		{
			options.kernelx               = (double *)mxMalloc(npt_ker*sizeof(double));
			for(i = 0 ; i < npt_ker ; i++)
			{
				options.kernelx[i]        = kernelx_default[i];
			}	
		}

		mxtemp                            = mxGetField( prhs[1], 0, "kernely" );
		if(mxtemp != NULL)
		{
			options.kernely               = mxGetPr(mxtemp);
			options.kyy                   = (int) mxGetM(mxtemp);
			options.kxy                   = (int) mxGetN(mxtemp);
		}
		else
		{
			options.kernely               = (double *)mxMalloc(npt_ker*sizeof(double));
			for(i = 0 ; i < npt_ker ; i++)
			{
				options.kernely[i]        = kernely_default[i];
			}	
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "deltax");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			dummy                         = (int) floor(nx*scalemin);
			if( (tempint < 1) || (tempint > dummy) )
			{
				mexPrintf("deltax must be > 0 and floor(nx*min(scale)), force to 1");	
				options.deltax           = dummy;
			}
			else
			{
				options.deltax         = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "deltay");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			dummy                         = (int) floor(ny*scalemin);
			if( (tempint < 1))
			{
				mexPrintf("deltay must be > 0 and floor(ny*min(scale)), force to 1");	
				options.deltay           = dummy;
			}
			else
			{
				options.deltay         = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "n");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 2) )
			{
				mexPrintf("n > 1, force to 2\n");	
				options.n                = 2;
			}
			else
			{
				options.n                = tempint;
			}
		}

		mxtemp                            = mxGetField( prhs[1], 0, "ldporient" );
		if(mxtemp != NULL)
		{
			options.norient               = (int) mxGetN(mxtemp);
			options.ldporient             = (int *)mxMalloc(npt_ldp*sizeof(int));
			tmp                           = mxGetPr(mxtemp);
			for(i = 0 ; i < npt_ldp ; i++)
			{
				options.ldporient[i]      = (int)tmp[i];
			}	
		}
		else
		{
			options.ldporient               = (int *)mxMalloc(npt_ldp*sizeof(int));
			for(i = 0 ; i < npt_ldp ; i++)
			{
				options.ldporient[i]        = ldporient_default[i];
			}	
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "color");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 0) || (tempint > 6))
			{
				mexPrintf("color = {0,1,2,3,4,5}, force to 0\n");	
				options.color             = 0;
			}
			else
			{
				options.color            = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "patchsize");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 1))
			{
				mexPrintf("patchsize > 0, force to 9");	
				options.patchsize         = 9;
			}
			else
			{
				options.patchsize         = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "improvedLDP");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 0) || (tempint > 1))
			{
				mexPrintf("improvedLDP = {0,1}, force to 0\n");	
				options.improvedLDP      = 0;
			}
			else
			{
				options.improvedLDP      = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "rmextremebins");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 0) || (tempint > 1))
			{
				mexPrintf("rmextremebins = {0,1}, force to 0\n");	
				options.rmextremebins     = 0;
			}
			else
			{
				options.rmextremebins     = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "norm");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			tempint                       = (int) tmp[0];
			if( (tempint < 0) || (tempint > 4) )
			{
				mexPrintf("norm = {0,1,2,3,4}, force to 0\n");	
				options.norm             = 0;
			}
			else
			{
				options.norm             = tempint;
			}
		}

		mxtemp                            = mxGetField(prhs[1] , 0 , "clamp");
		if(mxtemp != NULL)
		{
			tmp                           = mxGetPr(mxtemp);
			temp                          = tmp[0];

			if( (temp < 0.0) )
			{
				mexPrintf("clamp must be >= 0, force to 0.2\n");	
				options.clamp             = 0.2;
			}
			else
			{
				options.clamp             = temp;
			}
		}
	}
	else
	{
		options.nscale                = 1;
		options.scale                 = (double *)mxMalloc(1*sizeof(double));
		options.scale[0]              = 1.0;

		options.kernelx               = (double *)mxMalloc(npt_ker*sizeof(double));
		options.kernely               = (double *)mxMalloc(npt_ker*sizeof(double));

		for(i = 0 ; i < npt_ker ; i++)
		{
			options.kernelx[i]        = kernelx_default[i];
			options.kernely[i]        = kernely_default[i];
		}	
	}

	if((mxGetNumberOfDimensions(prhs[0]) == 2))
	{
		options.color  = 0;
		I              = (double *)mxMalloc(nynx*sizeof(double));
		dimcolor       = 1;
		for (i = 0 ; i < nynx ; i++)
		{
			I[i]       = (double)im[i];
		}
	}
	else
	{
		if((options.color == 0) )
		{
			I        = (double *)mxMalloc(nynx*sizeof(double));
			rgb2gray(im , ny , nx , I);
			dimcolor = 1;
		}
		else if (options.color == 1)
		{
			I        = (double *)mxMalloc(3*nynx*sizeof(double));
			for (i = 0 ; i < 3*nynx ; i++)
			{
				I[i] = (double)im[i];
			}
			dimcolor = 3;
		}
		else if (options.color == 2)
		{
			I        = (double *)mxMalloc(3*nynx*sizeof(double));
			rgb2nrgb(im , ny , nx , I);
			dimcolor = 3;
		}
		else if (options.color == 3)
		{
			I        = (double *)mxMalloc(3*nynx*sizeof(double));
			rgb2opponent(im , ny , nx , I);
			dimcolor = 3;
		}
		else if(options.color == 4)
		{
			I        = (double *)mxMalloc(2*nynx*sizeof(double));
			rgb2nopponent(im , ny , nx , I);
			dimcolor = 2;
		}
		else if(options.color == 5)
		{
			I        = (double *)mxMalloc(nynx*sizeof(double));
			rgb2hue(im , ny , nx , I);
			dimcolor = 1;
		}
	}

	/*----------------------- Outputs -------------------------------*/

	npts               =  options.deltay*options.deltax*options.nscale*dimcolor;

	plhs[0]            =  mxCreateDoubleMatrix(Nbins*options.norient*(1+options.improvedLDP) , npts , mxREAL);
	feat               =  mxGetPr(plhs[0]);

	plhs[1]            =  mxCreateDoubleMatrix(6 , npts , mxREAL);
	des                =  mxGetPr(plhs[1]);

	/*------------------------ Main Call ----------------------------*/

	denseMBLDP(I ,  ny , nx  , options , feat , des , npts , dimcolor );

	/*--------------------------- Free memory -----------------------*/

	if ( (nrhs > 1) && !mxIsEmpty(prhs[1]) )
	{
		if ( (mxGetField( prhs[1] , 0 , "scale" )) == NULL )
		{
			mxFree(options.scale);
		}
		if ( (mxGetField( prhs[1] , 0 , "kernelx" )) == NULL )
		{
			mxFree(options.kernelx);
		}
		if ( (mxGetField( prhs[1] , 0 , "kernely" )) == NULL )
		{
			mxFree(options.kernely);
		}
	}
	else
	{
		mxFree(options.scale);
		mxFree(options.kernelx);
		mxFree(options.kernely);
	}
	mxFree(options.ldporient);
	mxFree(I);
}

#else


#endif

/*----------------------------------------------------------------------------------------------------------------------------------------- */
void denseMBLDP(double *I , int ny , int nx , struct opts options , double *feat , double *des  , int npts , int dimcolor)
{
	double *kernelx = options.kernelx , *kernely = options.kernely;
	double *scale = options.scale;
	double *grady, *gradx , *tempgrady , *tempgradx , *mag;
	double *II , *IItemp;
	int i , l , m , v , s;
	int nynx = ny*nx , vnynx , coscale = 0 , nscale = options.nscale , improvedLDP = options.improvedLDP;
	int kyy = options.kyy , kxy = options.kxy , kyx = options.kyx , kxx = options.kxx ;
	int patchsize = options.patchsize , half_patch = (patchsize - 1)/2 , half_patch1 = -half_patch + 1, scalecurrent;
    int norient = options.norient , NBINS = 256 ,  Nbins = NBINS*norient;
	int deltay=options.deltay, deltax=options.deltax;
	double tempx , tempy , temp;
	int origy, origx;
	int nygrady , nxgrady , nygradx , nxgradx , noconv ;
	int nynew , nxnew  , lynew , lxnew , offsety , offsetx , co  , indd;

	if(improvedLDP)
	{
		Nbins          *= 2;
	}

	mag         = (double *)malloc(nynx*sizeof(double));
	II          = (double *)malloc(nynx*sizeof(double));
	IItemp      = (double *)malloc(nynx*sizeof(double));

	if(((kyy == 1) && (kxy == 1) && (kyx == 1) && (kxx == 1)) )
	{
		noconv      = 1;
	}
	else
	{
		noconv      = 0;
		nygrady     = ny + kyy - 1;
		nxgrady     = nx + kxy - 1;
		nygradx     = ny + kyx - 1;
		nxgradx     = nx + kxx - 1;

		grady       = (double *)malloc(nynx*sizeof(double));
		gradx       = (double *)malloc(nynx*sizeof(double));
		tempgrady   = (double *)malloc(nygrady*nxgrady*sizeof(double));
		tempgradx   = (double *)malloc(nygradx*nxgradx*sizeof(double));
	}

	co          = 0;
	for (v = 0 ; v < dimcolor ; v++)
	{
		vnynx    = v*nynx;
		if(noconv )
		{
			for (i = 0 ; i < nynx ; i++)
			{
				mag[i] = I[i + vnynx];
			}
		}
		else
		{
			for(i = 0 ; i < nygrady*nxgrady ; i++)
			{
				tempgrady[i] = 0.0;
			}
			conv2(I + vnynx , kernely , ny , nx , kyy , kxy , nygrady , nxgrady , grady , tempgrady);
			for(i = 0 ; i < nygradx*nxgradx ; i++)
			{
				tempgradx[i] = 0.0;
			}
			conv2(I + vnynx , kernelx , ny , nx , kyx , kxx , nygradx , nxgradx , gradx , tempgradx);

			for(i = 0 ; i < nynx ; i++)
			{
				tempy  = grady[i];
				tempx  = gradx[i];
				temp   = sqrt(tempy*tempy + tempx*tempx);
				mag[i] = temp;
			}
		}

		MakeIntegralImage(mag , II , nx, ny , IItemp);

		for (s = 0 ; s < nscale ; s++)
		{
			scalecurrent = (int) scale[s];
			nynew        = max(1 , ny - patchsize);
			nxnew        = max(1 , nx - patchsize);
			lynew        = (int) floor(nynew/(deltay - 1) );
			lxnew        = (int) floor(nxnew/(deltax - 1) );
			offsety      = max(0 , (int) ( floor((nynew - (deltay - 1)*lynew)/2))) + half_patch - 1;
			offsetx      = max(0 , (int) ( floor((nxnew - (deltax - 1)*lxnew)/2))) + half_patch - 1;

			for(l = 0 ; l < deltax ; l++) /* Loop shift on x-axis  */
			{
				origx  = offsetx + l*lxnew;
				for(m = 0 ; m < deltay ; m++)   /* Loop shift on y-axis  */
				{
					origy         = offsety + m*lynew;
					indd          = co*6;	

					compute_hmbldp(II , ny , nx , patchsize , patchsize , origy +  half_patch1 , origx +  half_patch1 , scalecurrent , feat + coscale , NBINS , options);

					des[0 + indd] = (double)(origy + 1);
					des[1 + indd] = (double)(origx + 1);
					des[2 + indd] = (double) scalecurrent;
					des[3 + indd] = (double) (v+1);
					des[4 + indd] = (double) (ny);
					des[5 + indd] = (double) (nx);
					co++;
					coscale      += Nbins ;
				}
			}
		}
	}

	free(mag);
	free(II);
	free(IItemp);
	if(!noconv )
	{
		free(tempgrady);
		free(tempgradx);
		free(grady);
		free(gradx);
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
void compute_hmbldp(double *II , int ny , int nx , int sy , int sx , int origy , int origx ,  int currentscale  , double *H , int NBINS , struct opts options )
{
	int i , a , xc , yc , xnw , ynw , xse , yse , sNbins;
	int norm = options.norm , n = options.n , n1 = n - 1;
	int improvedLDP = options.improvedLDP , rmextremebins = options.rmextremebins ;
	int norient = options.norient , *ldporient =options.ldporient , orient;
	int Nbins = NBINS;
	int startxc, startyc, endxc, endyc;
	double clamp = options.clamp;
	double sum, temp;
	double Ac , tmpA , sumA;
	unsigned short int valF;

	if(improvedLDP)
	{
		Nbins      = NBINS*2;
	}

	sNbins        = 0;

	for(a = 0 ; a < norient ; a++)
	{
		orient = ldporient[a];
		if(orient == 0)
		{
			startxc = 0;
			endxc   = n1;
			startyc = 0;
			endyc   = 0;
		}
		else if(orient == 1)
		{
			startxc = 0;
			endxc   = n1;
			startyc = 0;
			endyc   = n1;
		}
		else if(orient == 2)
		{
			startxc = 0;
			endxc   = 0;
			startyc = 0;
			endyc   = n1;
		}
		else if(orient == 3)
		{
			startxc = n1;
			endxc   = 0;
			startyc = 0;
			endyc   = n1;
		}

		for (xc = (1 + startxc)*currentscale + origx ; xc < sx - (2 + endxc)*currentscale + origx ; xc++)
		{
			for (yc = (1 + startyc)*currentscale + origy ; yc < sy - (2 + endyc)*currentscale + origy ; yc++)
			{
				xnw   = xc - currentscale;
				ynw   = yc - currentscale;

				xse   = xc + currentscale;
				yse   = yc + currentscale;

				Ac    = ld(II , xc  , yc  , currentscale , n1 , a , ny);

				valF  = 0;
				sumA  = Ac;

				tmpA  = ld(II , xnw , ynw , currentscale , n1 , a , ny);
				if(tmpA*Ac <= 0)
				{
					valF |= 0x01;
				}

				sumA += tmpA;

				tmpA  = ld(II , xc  , ynw , currentscale , n1 , a , ny);
				if(tmpA*Ac <= 0)
				{
					valF |= 0x02;
				}

				sumA += tmpA;

				tmpA  = ld(II , xse , ynw , currentscale , n1 , a , ny);
				if(tmpA*Ac <= 0)
				{
					valF |= 0x04;
				}

				sumA += tmpA;

				tmpA  = ld(II , xse , yc  , currentscale , n1 , a , ny);
				if(tmpA*Ac <= 0)
				{
					valF |= 0x08;
				}

				sumA += tmpA;

				tmpA  = ld(II , xse , yse , currentscale , n1 , a , ny);
				if(tmpA*Ac <= 0)
				{
					valF |= 0x10;
				}

				sumA += tmpA;

				tmpA  =  ld(II , xc  , yse , currentscale , n1 , a , ny);
				if(tmpA*Ac <= 0)
				{
					valF |= 0x20;
				}
				sumA += tmpA;

				tmpA =  ld(II , xnw , yse , currentscale , n1 , a , ny);
				if(tmpA*Ac <= 0)
				{
					valF |= 0x40;
				}
				sumA += tmpA;

				tmpA =  ld(II , xnw , yc , currentscale , n1 , a , ny);
				if(tmpA*Ac <= 0)
				{
					valF |= 0x80;
				}

				sumA += tmpA;

				if(improvedLDP)
				{
					if(sumA > 9*Ac)
					{
						valF |= 0x100; 
					}
				}
				H[valF + sNbins]++;
			}
		}


		if(rmextremebins == 1)
		{
			if(improvedLDP)
			{
				H[0 + sNbins] = H[Nbins-1 + sNbins] = H[Nbins + sNbins] = H[NBINS-1 + sNbins] = 0.0;
			}
			else
			{
				H[0 + sNbins] = H[Nbins-1 + sNbins] = 0.0;
			}
		}
			if(norm == 1)
			{
				sum       = 0.0;
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					sum   += H[i];
				}

				sum = 1.0/(sum + tiny);
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					H[i]   *= sum;
				}
			}
			if(norm == 2)
			{
				sum       = 0.0;
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					temp   = H[i];
					sum   += temp*temp;
				}

				sum = 1.0/sqrt(sum + verytiny);
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					H[i]   *= sum;
				}
			}
			if(norm == 3)
			{
				sum       = 0.0;
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					sum   += H[i];
				}

				sum = 1.0/(sum + tiny);
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					H[i]   = sqrt(H[i]*sum);
				}
			}
			if(norm == 4)
			{
				sum       = 0.0;
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					temp   = H[i];
					sum   += temp*temp;
				}

				sum = 1.0/sqrt(sum + verytiny);
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					H[i]   *= sum;
					if(H[i] > clamp)
					{
						H[i] = clamp;
					}
				}

				sum       = 0.0;
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					temp   = H[i];
					sum   += temp*temp;
				}

				sum = 1.0/sqrt(sum + verytiny);
				for(i = sNbins ; i < (NBINS + sNbins) ; i++)
				{
					H[i]   *= sum;
				}

			}		
		sNbins  += NBINS;
	}
}

/*----------------------------------------------------------------------------------------------------------------------------------------- */
void MakeIntegralImage(double *pIn, double *pOut, int iXmax, int iYmax , double *pTemp)
{
	/* Variable declaration */
	int x , y , indx = 0;

	for(x=0 ; x<iXmax ; x++)
	{
		pTemp[indx]     = pIn[indx];
		indx           += iYmax;
	}

	for(y = 1 ; y<iYmax ; y++)
	{
		pTemp[y]        = pTemp[y - 1] + pIn[y];
	}

	pOut[0]             = pIn[0];
	indx                = iYmax;

	for(x=1 ; x<iXmax ; x++)
	{
		pOut[indx]      = pOut[indx - iYmax] + pTemp[indx];
		indx           += iYmax;
	}
	for(y = 1 ; y<iYmax ; y++)
	{
		pOut[y]         = pOut[y - 1] + pIn[y];
	}

	/* Calculate integral image */

	indx                = iYmax;
	for(x = 1 ; x < iXmax ; x++)
	{
		for(y = 1 ; y < iYmax ; y++)
		{
			pTemp[y + indx]    = pTemp[y - 1 + indx] + pIn[y + indx];
			pOut[y + indx]     = pOut[y + indx - iYmax] + pTemp[y + indx];
		}
		indx += iYmax;
	}
}
/*----------------------------------------------------------------------------------------------------------------------------------------------*/
double ld(double *II , int x , int y , int s , int n , int a , int Ny)
{
	if (n > 1)
	{
		if(a == 0)
		{
			return (ld(II , x , y , s , n - 1 , a , Ny) - ld(II , x + s , y , s , n - 1 , a , Ny));
		}
		else if(a == 1)
		{
			return (ld(II , x , y , s , n - 1 , a , Ny) - ld(II , x + s , y + s , s , n - 1 , a , Ny));
		}
		else if(a == 2)
		{
			return (ld(II , x , y , s , n - 1 , a , Ny) - ld(II , x  , y + s , s , n - 1 , a , Ny));
		}
		else
		{
			return (ld(II , x , y , s , n - 1 , a , Ny) - ld(II , x - s , y + s , s , n - 1 , a , Ny));
		}
	}
	else
	{
		if(a == 0)
		{
			return (Area(II , x , y , s  , s  , Ny) - Area(II , x + s , y , s  , s  , Ny));
		}
		else if(a == 1)
		{
			return (Area(II , x , y , s  , s  , Ny) - Area(II , x + s , y + s , s  , s  , Ny));
		}
		else if(a == 2)
		{
			return (Area(II , x , y , s  , s  , Ny) - Area(II , x  , y + s , s  , s  , Ny));
		}
		else
		{
			return (Area(II , x , y , s  , s  , Ny) - Area(II , x - s , y + s , s  , s  , Ny));
		}
	}
}

/*----------------------------------------------------------------------------------------------------------------------------------------------*/
double Area(double *II , int x , int y , int w , int h , int Ny)
{
	int h1 = h-1 , w1 = w-1 , x1 = x-1, y1 = y-1;

	if( (x == 0) && (y==0))
	{
		return (II[h1 + w1*Ny]);
	}
	if( (x==0) ) 
	{
		return(II[(y+h1) + w1*Ny] - II[y1 + w1*Ny]);
	}
	if( (y==0) ) 
	{
		return(II[h1 + (x+w1)*Ny] - II[h1 + x1*Ny]);
	}
	else
	{
		return (II[(y+h1) + (x+w1)*Ny] - (II[y1 + (x+w1)*Ny] + II[(y+h1) + x1*Ny]) + II[y1 + x1*Ny]);
	}
}
/*----------------------------------------------------------------------------------------------------------------------------------------- */
void conv2(double *A, double *B, int nyA , int nxA , int nyB , int nxB , int nyD , int nxD , double *C , double *D)
{
	int p , q , index; 
	int lengthA = nyA*nxA , lengthB = nyB*nxB , lengthD = nyD*nxD , nyA1 = nyA - 1 , nyB1 = nyB - 1;
	int nyAtemp , nxAtemp ,  nyBtemp , nxBtemp , nyDtemp , nxDtemp;
	int nytemp = (int)ceil((1.0 + nyB)*0.5) - 1 , nxtemp = (int)ceil((1.0 + nxB)*0.5) - 1;

	for( p = 0 ; p < lengthD ; p++)
	{
		D[p] = 0.0;
	}

	nyAtemp  = -1;
	nxAtemp  = 0;

	for (p = 0 ; p < lengthA ; p++) 
	{
		nyAtemp  += 1;

		if (nyAtemp > nyA1)
		{ 
			nyAtemp    = 0; 
			nxAtemp   += 1; 
		}

		nyBtemp  = -1;
		nxBtemp  = 0;
		for (q = 0 ; q < lengthB ; q++) 
		{
			nyBtemp  += 1;

			if (nyBtemp > nyB1)
			{ 
				nyBtemp    = 0; 
				nxBtemp   += 1; 				
			}

			index           = nyAtemp + nyBtemp + (nxAtemp + nxBtemp)*nyD;
			D[index]       += A[p]*B[q];
		}
	}	

	nyDtemp = nytemp - 1;	 
	nxDtemp = nxtemp;
	nyA1   += nytemp;

	for (p = 0 ; p < lengthA ; p++) 
	{
		nyDtemp  += 1;

		if(nyDtemp > nyA1)
		{
			nyDtemp  = nytemp;
			nxDtemp += 1;
		}

		index     = nyDtemp + nxDtemp*nyD;
		C[p]      = D[index];
	}
}

/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2gray(unsigned char *rgb , int Ny , int Nx , double *gray)
{
	/* T    = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);

	vect = T(1 , :);

	*/

	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx;

	for (i = 0 ; i  < NyNx ; i++)
	{
		gray[i]   = 0.298936021293776*rgb[i] + 0.587043074451121*rgb[i + NyNx] + 0.114020904255103*rgb[i + NyNx2];
	}
}

/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2nrgb(unsigned char *rgb , int Ny , int Nx , double *nrgb)
{
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx ;
	double r , g , b , cte;

	for (i = 0 ; i  < NyNx ; i++)
	{
		r                = (double)rgb[i];
		g                = (double)rgb[i + NyNx];
		b                = (double)rgb[i + NyNx2];

		cte              = (r + g + b);
		if(cte != 0.0)
		{
			cte          = 1.0/cte;
		}
		else
		{
			cte          = 1.0;
		}

		nrgb[i]          = r*cte;
		nrgb[i + NyNx]   = g*cte;
		nrgb[i + NyNx2]  = b*cte;
	}
}

/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2opponent(unsigned char *rgb , int Ny , int Nx , double *opponent)
{
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx ;
	double r , g , b ;

	for (i = 0 ; i  < NyNx ; i++)
	{
		r                    = (double)rgb[i];
		g                    = (double)rgb[i + NyNx];
		b                    = (double)rgb[i + NyNx2];

		opponent[i]          = (r - g)*invsqrt2;
		opponent[i + NyNx]   = (r + g - 2*b)*invsqrt6;
		opponent[i + NyNx2]  = (r + g + b)*invsqrt3;
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2nopponent(unsigned char *rgb , int Ny , int Nx , double *nopponent)
{
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx;
	double r , g , b , cte ;

	for (i = 0 ; i  < NyNx ; i++)
	{
		r                     = (double)rgb[i];
		g                     = (double)rgb[i + NyNx];
		b                     = (double)rgb[i + NyNx2];

		cte                   = (r + g + b);
		if(cte != 0.0)
		{
			cte               = 1.0/(cte*invsqrt2);
		}
		else
		{
			cte               = 1.0;
		}

		nopponent[i]          = sqrt3*(r - g)*cte;
		nopponent[i + NyNx]   = (r + g - 2*b)*cte;
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
void rgb2hue(unsigned char *rgb , int Ny , int Nx , double *hue)
{
	int i , NyNx = Ny*Nx , NyNx2 = 2*NyNx ;
	double r , g , b ;

	for (i = 0 ; i  < NyNx ; i++)
	{
		r         = (double)rgb[i];
		g         = (double)rgb[i + NyNx];
		b         = (double)rgb[i + NyNx2];
		hue[i]    = atan(sqrt3*(r - g)/(r + g - 2*b + verytiny));
	}
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/

