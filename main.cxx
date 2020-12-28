#include <iostream>
#include <string>

#include "json/json.h"
#include <fstream> 
#include <vtkMetaImageReader.h>
#include <vtkPlane.h>

#include "itkImage.h"
#include "itkVTKImageToImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkRayCastInterpolateImageFunction.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileWriter.h"

const unsigned int Dimension = 3;
//input image
typedef float PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
//output image
typedef unsigned char PngPixelType;
typedef itk::Image< PngPixelType, 2 > OutputImageType;

struct feature_3d_t
{
	float x;
	float y;
	float z;
	std::string name;
};

struct feature_2d_t
{
	int i;
	int j;
	std::string name;
	std::string DRR_filename;
};
///////////////////////////////////////////////////////////////////////////////////////////////////
void generate_DRR(const float rx, const float ry, const float rz,
				  const float tx, const float ty, const float tz,
				  const float d_f3d,
				  const float d_3d2d,
			  	  const int npixel_x, const int npixel_y,
				  const float res_x, const float res_y,
				  const ImageType::Pointer input_image,
				  const std::vector<feature_3d_t> features_3d,
				  const std::string output_filename,
				  std::vector<feature_2d_t> &features_2d)
{
	std::cout << "generate_DRR() - begin" << std::endl;
	// Don't move drr, move head!!!

	float cx = 0.0;
	float cy = 0.0;
	float cz = 0.0;
	double threshold = 0.0;
	float o2Dx = 0.0;
	float o2Dy = 0.0;

    /* Creation of a ResampleImageFilter enables coordinates for
	   each of the pixels in the DRR image to be generated. These
	   coordinates are used by the RayCastInterpolateImageFunction
	   to determine the equation of each corresponding ray which is cast
	   through the input volume.
	*/
	typedef itk::ResampleImageFilter<ImageType, ImageType > DRRFilterType;
	DRRFilterType::Pointer drr_filter = DRRFilterType::New();
	drr_filter->SetDefaultPixelValue(0);
	
    /* An Euler transformation is defined to position the input volume.
	   The ResampleImageFilter uses this transform to position the
	   output DRR image for the desired view.
	*/
	typedef itk::CenteredEuler3DTransform< double >  TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->SetComputeZYX(true);
	TransformType::OutputVectorType translation;
	translation[0] = 0.0;
    translation[1] = 0.0;
    translation[2] = 0.0;
   
	// Translation and Rotation matrix
	transform->SetTranslation(translation);
	transform->SetRotation(0.0, 0.0, 0.0);
	
	//Here is the origin of the lower left corner of the image
	ImageType::PointType   imOrigin = input_image->GetOrigin();
	std::cout << "imOrigin: " << imOrigin[0] << "," << imOrigin[1] << "," << imOrigin[2] << std::endl;
	//Here is the real physical pixel width of the image (mm)
	ImageType::SpacingType imRes = input_image->GetSpacing();
	std::cout << "imRes: " << imRes[0] << "," << imRes[1] << "," << imRes[2] << std::endl;
	
	typedef ImageType::RegionType     ImageRegionType;
	typedef ImageRegionType::SizeType ImageSizeType;
	ImageRegionType imRegion = input_image->GetBufferedRegion();
	//Get image size
	ImageSizeType   imSize = imRegion.GetSize();
	std::cout << "imSize: " << imSize[0] << "," << imSize[1] << "," << imSize[2] << std::endl;
	//Calculate the physical coordinates of the center point of 
	//the image with the lower left corner as the origin
	imOrigin[0] += imRes[0] * static_cast<double>(imSize[0]) / 2.0;
	imOrigin[1] += imRes[1] * static_cast<double>(imSize[1]) / 2.0;
	imOrigin[2] += imRes[2] * static_cast<double>(imSize[2]) / 2.0;
	std::cout << "imOrigin - updated: " << imOrigin[0] << "," << imOrigin[1] << "," << imOrigin[2] << std::endl;
	//The total center coordinates of the image, if not set, it is cx=0,cy=0,cz=0
	TransformType::InputPointType center;
	center[0] = cx + imOrigin[0];
	center[1] = cy + imOrigin[1];
	center[2] = cz + imOrigin[2];
	//Set the physical coordinates of the conversion center
	transform->SetCenter(center);
	std::cout << "center: " << center[0] << "," << center[1] << "," << center[2] << std::endl;

    /* The RayCastInterpolateImageFunction is instantiated and passed the transform
	   object. The RayCastInterpolateImageFunction uses this transform to reposition 
	   the x-ray source such that the DRR imageand x-ray source move as one around
	   the input volume. This coupling mimics the rigid geometry of the x-ray gantry.
	*/
	typedef itk::RayCastInterpolateImageFunction<ImageType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetTransform(transform);

	// specify a threshold above which the volume's intensities will be integrated.
	interpolator->SetThreshold(threshold);

	/* The ray-cast interpolator needs to know the initial position of the
	   ray source or focal point. In this example we place the input
	   volume at the origin and a position between the ray source and the
	   screen. The distance between the ray source and the 3d object is 
	   the "source to 3d image distance" d_f3d and is specified by the user.
    */
	// set x-ray source
	InterpolatorType::InputPointType focalpoint;
	focalpoint[0] = imOrigin[0];
	focalpoint[1] = imOrigin[1];
	focalpoint[2] = imOrigin[2] - d_f3d;
	interpolator->SetFocalPoint(focalpoint);
	std::cout << "focalpoint: " << focalpoint[0] << "," << focalpoint[1] << "," << focalpoint[2] << std::endl;
	
    // Having initialised the interpolator we pass the object to the resample filter.
	drr_filter->SetInterpolator(interpolator);
	drr_filter->SetTransform(transform);

	//Transform original image
	std::vector<ImageType::PointType> features_xray_3d;
	{
		typedef itk::ResampleImageFilter<ImageType, ImageType > TransformFilterType;
		TransformFilterType::Pointer transformation_filter = TransformFilterType::New();
		
		TransformType::Pointer transform2 = TransformType::New();
		transform2->SetComputeZYX(true);
		TransformType::OutputVectorType translation2;
		translation2[0] = tx;
		translation2[1] = ty;
		translation2[2] = tz;
		
        // constant for converting degrees into radians
        const double dtr = ( std::atan(1.0) * 4.0 ) / 180.0;
        
		transform2->SetTranslation(translation2);
		transform2->SetRotation(dtr*rx, dtr*ry, dtr*rz);

		ImageType::PointType   imOrigin2 = input_image->GetOrigin();
		ImageType::SpacingType imRes2 = input_image->GetSpacing();
		
		ImageRegionType imRegion2 = input_image->GetBufferedRegion();
		ImageSizeType   imSize2 = imRegion2.GetSize();
		
		imOrigin2[0] += imRes2[0] * static_cast<double>(imSize2[0]) / 2.0;
		imOrigin2[1] += imRes2[1] * static_cast<double>(imSize2[1]) / 2.0;
		imOrigin2[2] += imRes2[2] * static_cast<double>(imSize2[2]) / 2.0;
		
		TransformType::InputPointType center2;
		center2[0] = cx + imOrigin2[0];
		center2[1] = cy + imOrigin2[1];
		center2[2] = cz + imOrigin2[2];
		std::cout << "center2: " << center2[0] << "," << center2[1] << "," << center2[2] << std::endl;
		transform2->SetCenter(center2);
		cout << "transform2 set" << endl;

		TransformType::InverseTransformBasePointer transform2_inv = transform2->GetInverseTransform();

		typedef itk::LinearInterpolateImageFunction< ImageType, double >  LinearInterpolatorType;
		LinearInterpolatorType::Pointer interpolator2 = LinearInterpolatorType::New();
		transformation_filter->SetInterpolator(interpolator2);
		transformation_filter->SetDefaultPixelValue(-1024);
		transformation_filter->SetInput(input_image);
		transformation_filter->SetTransform(transform2);

		const ImageType::SpacingType & spacing = input_image->GetSpacing();
		const ImageType::PointType & origin = input_image->GetOrigin();
		ImageType::SizeType size = input_image->GetLargestPossibleRegion().GetSize();
		transformation_filter->SetOutputOrigin(origin);
		transformation_filter->SetOutputSpacing(spacing);
		transformation_filter->SetOutputDirection(input_image->GetDirection());
		transformation_filter->SetSize(size);
		transformation_filter->Update();
		cout << "filter updated" << endl;

		drr_filter->SetInput(transformation_filter->GetOutput());
		
		//Given a line defined by the two points p1,p2
		//and a plane defined by the normal n and point p0, compute an intersection.
		double p0[Dimension];
		p0[0] = tx + imOrigin[0] + o2Dx - res_x * (static_cast<double>(npixel_x) - 1.) / 2.;
		p0[1] = ty + imOrigin[1] + o2Dy - res_y * (static_cast<double>(npixel_y) - 1.) / 2.;
		p0[2] = tz + imOrigin[2] + d_3d2d;
		std::cout << "tx: " << tx << std::endl;
		std::cout << "ty: " << ty << std::endl;
		std::cout << "tz: " << tz << std::endl;
		std::cout << "p0: " << p0[0] << ", " << p0[1] << ", " << p0[2] << std::endl;
		// normal n
		double normal[Dimension];
		normal[0] = 0.0;
		normal[1] = 0.0;
		normal[2] = 1.0;

		double p1[Dimension];
		p1[0] = focalpoint[0];
		p1[1] = focalpoint[1];
		p1[2] = focalpoint[2];

		for (unsigned i = 0; i < features_3d.size(); ++i)
		{
			TransformType::InputPointType feature_pt;
			feature_pt[0] = features_3d[i].x;
			feature_pt[1] = features_3d[i].y;
			feature_pt[2] = features_3d[i].z;

			TransformType::InputPointType transformed_feature_pt;
			transformed_feature_pt = transform2_inv->TransformPoint(feature_pt);

			{
				double p2[Dimension];
				p2[0] = transformed_feature_pt[0];
				p2[1] = transformed_feature_pt[1];
				p2[2] = transformed_feature_pt[2];

				double t;//The parametric coordinate along the line 
				
				double feature_xray[Dimension];//the coordinates of intersection
				//A zero is returned if the plane and line do not intersect between (0<=t<=1). 
				//If the plane and line are parallel, zero is returned and t is set to VTK_LARGE_DOUBLE.
				vtkPlane::IntersectWithLine(p1, p2, normal, p0, t, feature_xray);
				std::cout << "xray_feature - 3d: " << feature_xray[0] << ", " << feature_xray[1] << ", " << feature_xray[2] << std::endl;
				std::cout << "t: " << t << std::endl;
				features_xray_3d.push_back(feature_xray);
			}
		}
	}

	// The size and resolution of the output DRR image is specified via the resample filter.
	ImageType::SizeType   size;//Number of pixels in X and Y directions
	size[0] = npixel_x;  // number of pixels along X of the 2D DRR image
	size[1] = npixel_y;  // number of pixels along Y of the 2D DRR image
	size[2] = 1;         // only one slice
	std::cout << "npixel_y: " << npixel_y << std::endl;
	drr_filter->SetSize(size);
	
	ImageType::SpacingType spacing;//Set the physical size of pixels
	spacing[0] = res_x;  // pixel spacing along X of the 2D DRR image [mm]
	spacing[1] = res_y;  // pixel spacing along Y of the 2D DRR image [mm]
	spacing[2] = 1.0;    // slice thickness of the 2D DRR image [mm]
	drr_filter->SetOutputSpacing(spacing);
	
    /* In addition the position of the DRR is specified. The default
	   position of the input volume, prior to its transformation is
	   half-way between the ray source and screen and unless specified
	   otherwise the normal from the "screen" to the ray source passes
	   directly through the centre of the DRR.
	*/
	double origin[Dimension];
	origin[0] = imOrigin[0] + o2Dx - res_x * (static_cast<double>(npixel_x) - 1.) / 2.;
	origin[1] = imOrigin[1] + o2Dy - res_y * (static_cast<double>(npixel_y) - 1.) / 2.;
	origin[2] = imOrigin[2] + d_3d2d;
	std::cout << "origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;
	
	drr_filter->SetOutputOrigin(origin);
	drr_filter->Print(std::cout);
	drr_filter->Update();

	cout << "drr_filter updated" << endl;
	

	// Save DRR as png
	std::string drr_filename = "PNG/"+output_filename + ".png";
	{
		typedef itk::RescaleIntensityImageFilter<ImageType, ImageType > RescaleFilterType;
		RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
		rescaler->SetOutputMinimum(0);
		rescaler->SetOutputMaximum(255);
		rescaler->SetInput(drr_filter->GetOutput());

		typedef itk::ExtractImageFilter< ImageType, OutputImageType > SlicerType;
		SlicerType::Pointer slicer = SlicerType::New();
		ImageType::RegionType inputRegion = drr_filter->GetOutput()->GetLargestPossibleRegion();
		ImageType::SizeType size = inputRegion.GetSize();
		size[2] = 0;
		cout << "size - slicer: " << size[0] << ", " << size[1] << ", " << size[2] << endl;


		ImageType::IndexType start;
		const unsigned int sliceNumber = 0;
		start[0] = 0;
		start[1] = 0;
		start[2] = sliceNumber;
		ImageType::RegionType desiredRegion;
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(start);
		cout << "index - slicer: " << start[0] << ", " << start[1] << ", " << start[2] << endl;

		slicer->SetExtractionRegion(desiredRegion);
		slicer->SetInput(rescaler->GetOutput());
		slicer->SetDirectionCollapseToSubmatrix();
		slicer->Update();
		cout << "slicer updated" << endl;

		typedef itk::ImageFileWriter< OutputImageType >  Writer2DType;
		Writer2DType::Pointer writer2d = Writer2DType::New();
		writer2d->SetInput(slicer->GetOutput());
		writer2d->SetFileName(drr_filename);

		try
		{
			std::cout << "Writing image: " << drr_filename << std::endl;
			writer2d->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << "ERROR: ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}
	}

	// Convert feature 3d coordinate to 2d pixel coordinates
	std::vector<OutputImageType::PointType> features_indices;
	for (unsigned i = 0; i< features_xray_3d.size(); ++i)
	{
		ImageType::PointType point;
		point[0] = features_xray_3d[i][0];
		point[1] = features_xray_3d[i][1];
		point[2] = features_xray_3d[i][2];

		ImageType::IndexType index_feature_3d;
		drr_filter->GetOutput()->TransformPhysicalPointToIndex(point, index_feature_3d);

		//Drop one dimension. The output of drr_filter is actually a one voxel thin image.
		double index_feature_2d[2];
		index_feature_2d[0] = index_feature_3d[0];
		index_feature_2d[1] = index_feature_3d[1];

		cout << "index_xray_feature - index: " << index_feature_3d[0] << ", " <<
			index_feature_3d[1] << ", " <<
			index_feature_3d[2] << endl;

		features_indices.push_back(index_feature_2d);
		feature_2d_t feature_2d;
		feature_2d.i = index_feature_2d[0];
		feature_2d.j = index_feature_2d[1];
		feature_2d.name = features_3d[i].name;
		feature_2d.DRR_filename = drr_filename;

		features_2d.push_back(feature_2d);
	}

	std::cout << "generate_DRR() - end" << std::endl;


}
///////////////////////////////////////////////////////////////////////////////////////////////////
void parse_features(const std::string filename,
	std::vector<feature_3d_t> &features)
{
	std::cout << "parse_features - begin" << std::endl;
	std::cout << "features size: " << features.size() << std::endl;
	std::ifstream file(filename);

	Json::Value root;
	file >> root;

	for (Json::Value::ArrayIndex i = 0; i != root["features"].size(); i++)
	{
		feature_3d_t feature;
		feature.name = root["features"][i]["name"].asString();
		feature.x = root["features"][i]["x"].asFloat();
		feature.y = root["features"][i]["y"].asFloat();
		feature.z = root["features"][i]["z"].asFloat();
		features.push_back(feature);
	}

	std::cout << "features size: " << features.size() << std::endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////
void parse_config(const std::string filename,
	float &tx_min, float &tx_max, float &tx_delta,
	float &ty_min, float &ty_max, float &ty_delta,
	float &tz_min, float &tz_max, float &tz_delta,
	float &rx_min, float &rx_max, float &rx_delta,
	float &ry_min, float &ry_max, float &ry_delta,
	float &rz_min, float &rz_max, float &rz_delta,
	float &d_f3d,float &d_3d2d, int &npixel_x, int &npixel_y, float &res_x, float &res_y
)
{
	std::ifstream file(filename);

	Json::Value root;
	file >> root;

	// Parse xray settings
	d_f3d = root["xray"]["d_f3d"].asFloat();
	d_3d2d = root["xray"]["d_3d2d"].asFloat();
	npixel_x = root["xray"]["npixel_x"].asInt();
	npixel_y = root["xray"]["npixel_y"].asInt();
	res_x = root["xray"]["res_x"].asFloat();
	res_y = root["xray"]["res_y"].asFloat();
	std::cout << "parsing - d_f3d: " << d_f3d << std::endl;
	std::cout << "parsing - d_3d2d: " << d_3d2d << std::endl;
	std::cout << "parsing - npixel_x: " << npixel_x << std::endl;
	std::cout << "parsing - npixel_y: " << npixel_y << std::endl;
	std::cout << "parsing - res_x: " << res_x << std::endl;
	std::cout << "parsing - res_y: " << res_y << std::endl;

	// Parse transformation
	tx_min = root["transformation"]["tx_min"].asFloat();
	tx_max = root["transformation"]["tx_max"].asFloat();
	tx_delta = root["transformation"]["tx_delta"].asFloat();
	ty_min = root["transformation"]["ty_min"].asFloat();
	ty_max = root["transformation"]["ty_max"].asFloat();
	ty_delta = root["transformation"]["ty_delta"].asFloat();
	tz_min = root["transformation"]["tz_min"].asFloat();
	tz_max = root["transformation"]["tz_max"].asFloat();
	tz_delta = root["transformation"]["tz_delta"].asFloat();

	rx_min = root["transformation"]["rx_min"].asFloat();
	rx_max = root["transformation"]["rx_max"].asFloat();
	rx_delta = root["transformation"]["rx_delta"].asFloat();
	ry_min = root["transformation"]["ry_min"].asFloat();
	ry_max = root["transformation"]["ry_max"].asFloat();
	ry_delta = root["transformation"]["ry_delta"].asFloat();
	rz_min = root["transformation"]["rz_min"].asFloat();
	rz_max = root["transformation"]["rz_max"].asFloat();
	rz_delta = root["transformation"]["rz_delta"].asFloat();

	std::cout << "parsing - tx_min: " << tx_min << std::endl;
	std::cout << "parsing - tx_max: " << tx_max << std::endl;
	std::cout << "parsing - tx_delta: " << tx_delta << std::endl;
	std::cout << "parsing - ty_min: " << ty_min << std::endl;
	std::cout << "parsing - ty_max: " << ty_max << std::endl;
	std::cout << "parsing - ty_delta: " << ty_delta << std::endl;
	std::cout << "parsing - tz_min: " << tz_min << std::endl;
	std::cout << "parsing - tz_max: " << tz_max << std::endl;
	std::cout << "parsing - tz_delta: " << tz_delta << std::endl;

	std::cout << "parsing - rx_min: " << rx_min << std::endl;
	std::cout << "parsing - rx_max: " << rx_max << std::endl;
	std::cout << "parsing - rx_delta: " << rx_delta << std::endl;
	std::cout << "parsing - ry_min: " << ry_min << std::endl;
	std::cout << "parsing - ry_max: " << ry_max << std::endl;
	std::cout << "parsing - ry_delta: " << ry_delta << std::endl;
	std::cout << "parsing - rz_min: " << rz_min << std::endl;
	std::cout << "parsing - rz_max: " << rz_max << std::endl;
	std::cout << "parsing - rz_delta: " << rz_delta << std::endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
	if (argc < 5) 
	{ // We expect 4 arguments: the program name, the source path and the destination path
		std::cerr << "Usage: " << argv[0] 
			<< "  config_filename CT_filename features_filename stem" << std::endl;
		return 1;
	}

	const std::string json_filename = argv[1];
	float tx_min = 0.0;
	float tx_max = 0.0;
	float tx_delta = 0.0;
	float ty_min = 0.0;
	float ty_max = 0.0;
	float ty_delta = 0.0;
	float tz_min = 0.0;
	float tz_max = 0.0;
	float tz_delta = 0.0;

	float rx_min = 0.0;
	float rx_max = 0.0;
	float rx_delta = 0.0;
	float ry_min = 0.0;
	float ry_max = 0.0;
	float ry_delta = 0.0;
	float rz_min = 0.0;
	float rz_max = 0.0;
	float rz_delta = 0.0;
	float d_3d2d = 0.0;
	float d_f3d = 0.0;
	int npixel_x = 0;
	int npixel_y = 0;
	float res_x = 0;
	float res_y = 0;

	// Parse config file
	parse_config(json_filename, 
		tx_min, tx_max, tx_delta, ty_min, ty_max, ty_delta, tz_min, tz_max, tz_delta,
		rx_min, rx_max, rx_delta, ry_min, ry_max, ry_delta, rz_min, rz_max, rz_delta,
		d_f3d,d_3d2d, npixel_x, npixel_y, res_x, res_y
		);

	// Load CT image
	const std::string CT_filename = argv[2];
	
	vtkSmartPointer<vtkMetaImageReader> reader = vtkSmartPointer<vtkMetaImageReader>::New();
	reader->SetFileName(CT_filename.c_str());
	reader->Update();

	typedef itk::VTKImageToImageFilter< ImageType > VtkImporterType;
	VtkImporterType::Pointer vtk_importer = VtkImporterType::New();
	vtk_importer->SetInput(reader->GetOutput());	
	vtk_importer->Update();
	ImageType::Pointer input_image = vtk_importer->GetOutput();
	std::cout << "Image imported from vtk" << std::endl;

	// Load features 
	std::vector<feature_3d_t> features;
	const std::string features_filename = argv[3];
	parse_features(features_filename, features);
    
    // set a name for each trial 
	const std::string stem = argv[4];
	unsigned image_id = 0;
	
	float tx = tx_min;
	float ty = ty_min;
	float tz = tz_min;
	float rx = rx_min;
	float ry = ry_min;
	float rz = rz_min;

	std::vector<feature_2d_t> features_2d;
	features_2d.clear();
	
	// Save map between transform and DRR
	std::ofstream DRR_file;
	std::string DRR_filename = stem + "_DRRs.txt";

	DRR_file.open(DRR_filename);

	while (tx <= tx_max)
	{
		float ty = ty_min;
		while (ty <= ty_max)
		{
			float tz = tz_min;
			while (tz <= tz_max)
			{
				float rx = rx_min;
				while (rx <= rx_max)
				{
					float ry = ry_min;
					while (ry <= ry_max)
					{
						float rz = rz_min;
						while (rz <= rz_max)
						{
							std::cout << "tx: " << tx << std::endl;
							std::cout << "ty: " << ty << std::endl;
							std::cout << "tz: " << tz << std::endl;
							std::cout << "rx: " << rx << std::endl;
							std::cout << "ry: " << ry << std::endl;
							std::cout << "rz: " << rz << std::endl;
							std::cout << "image_id: " << image_id << std::endl;
							std::cout << "" << std::endl;

										
							std::ostringstream ss;
							ss << std::setw(7) << std::setfill('0') << image_id;
							std::string output_filename = stem + ss.str();
							DRR_file << CT_filename << " " << tx << " " << ty << " " << tz << " " << rx << " " << ry << " " << rz << " " << output_filename << std::endl;

							generate_DRR(rx, ry, rz, tx, ty, tz,
								d_f3d,d_3d2d, npixel_x, npixel_y, res_x, res_y, 
								input_image,
								features,
								output_filename,
								features_2d);
							
							++image_id;

							if (rz_delta < 0.0)
								rz = rz_max + 1.0;
							else if (rz_delta>0.0)
								rz += rz_delta;
							else
								break;
						}
						if (ry_delta < 0.0)
							ry = ry_max + 1.0;
						else if (ry_delta > 0.0)
							ry += ry_delta;
						else
							break;
					}
					if (rx_delta < 0.0)
						rx = rx_max + 1.0;
					else if (rx_delta>0.0)
						rx += rx_delta;
					else
						break;
				}
				if (tz_delta < 0.0)
					tz = tz_max + 1.0;
				else if (tz_delta > 0.0)
					tz += tz_delta;
				else 
					break;
			}
			if (ty_delta < 0.0)
				ty = ty_max + 1.0;
			else if (ty_delta > 0.0)
				ty += ty_delta;
			else
				break;
		}
		if (tx_delta < 0.0)
			tx = tx_max + 1.0;
		else if (tx_delta > 0.0)
			tx += tx_delta;
		else 
			break;
	}

	DRR_file.close();

	// Save features to txt file
	std::ofstream feature_file;
	std::string feature_filename = stem + "_features.txt";
	
	feature_file.open(feature_filename);
	
	for (unsigned i = 0; i < features_2d.size(); ++i)
	{
		std::cout << "Feature: " << features_2d[i].i << ", " << features_2d[i].j << ", " << features_2d[i].name << ", " << features_2d[i].DRR_filename << std::endl;
		feature_file << features_2d[i].name << " " << features_2d[i].i << " " << features_2d[i].j << " " << features_2d[i].DRR_filename << std::endl;
	}

	feature_file.close();
  return 0;
}
