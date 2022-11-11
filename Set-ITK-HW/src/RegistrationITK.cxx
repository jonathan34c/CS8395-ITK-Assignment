#include <iostream>

// include stuff here
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

//#include <itkImageRegistrationMethod.h>
#include <itkMultiResolutionImageRegistrationMethod.h>

#include <itkAffineTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkRegularStepGradientDescentOptimizer.h>


#include <itkResampleImageFilter.h>

#include <itkCommand.h>
#include <itkBSplineTransform.h>


#include <itkImageRegistrationMethodv4.h>
#include <itkMattesMutualInformationImageToImageMetricv4.h>
 
#include <itkTimeProbesCollectorBase.h>
#include <itkMemoryProbesCollectorBase.h>

#include <itkStatisticsImageFilter.h>

#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkVersorRigid3DTransform.h>

#include <itkAffineTransform.h>


// things we could/should do differently
// if we are inter-subject, A_T1 to B_T1, affine isn't good enough, should be deformable
// if we are intra-subject, A_T1 to A_T2, sum of square differences is not a good cross-modality metric, try correlation or mutual information

// declare our image type 
typedef itk::Image < int, 3 > myImageType ;

// todo: implement this class
class MyObserverCommandClass : public itk::Command 
{
public:
  itkNewMacro ( MyObserverCommandClass ) ;

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType ;
  //typedef OptimizerType * OptimizerPointerType ;
  void Execute ( itk::Object *caller, const itk::EventObject & event ) override
  {
    // std::cout << "HelloWorld" << std::endl ;
    OptimizerType::Pointer optimizer = dynamic_cast < OptimizerType * > ( caller ) ;
    unsigned int iterationNumber = optimizer->GetCurrentIteration () ;
    double currentMetricValue = optimizer->GetValue () ;
    std::cout << iterationNumber << " " << currentMetricValue << std::endl ;
    return ;
  }

  void Execute ( const itk::Object *caller, const itk::EventObject & event ) override
  {
    std::cout << "Hello world  -const" << std::endl ;
    return ;
  }
} ;


class MyRegistrationObserverCommandClass : public itk::Command 
{
public:
  itkNewMacro ( MyRegistrationObserverCommandClass ) ;

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType ;
  typedef itk::MultiResolutionImageRegistrationMethod < myImageType, myImageType > RegistrationMethod ;
  //typedef OptimizerType * OptimizerPointerType ;
  void Execute ( itk::Object *caller, const itk::EventObject & event ) override
  {
    RegistrationMethod::Pointer regMethod = dynamic_cast < RegistrationMethod * > ( caller ) ;
    int currentLevel = regMethod->GetCurrentLevel() ;
    std::cout << "Starting registration level: " << currentLevel << std::endl ;

    // optimizer is no longer the caller but we can get to it thru the regmethod
    // need the optimizer so we can change reg settings per level, e.g. step length or nIterations
    OptimizerType::Pointer optimizer = dynamic_cast < OptimizerType * > ( regMethod->GetModifiableOptimizer () ) ;

    if ( currentLevel == 0 )
      {
	optimizer->SetNumberOfIterations ( 60 ) ;
	std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
	std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
	optimizer->SetMinimumStepLength ( 0 ) ;
	optimizer->SetMaximumStepLength ( 0.125 ) ; // looks ok enough
      }
    else if ( currentLevel == 1 )
      {
	optimizer->SetNumberOfIterations ( 20 ) ;
	std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
	std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
	optimizer->SetMinimumStepLength ( 0 ) ;
	optimizer->SetMaximumStepLength ( 0.0625 ) ;
      }
    else
      {
	optimizer->SetNumberOfIterations ( 10 ) ;
	std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
	std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
	optimizer->SetMinimumStepLength ( 0 ) ;
	optimizer->SetMaximumStepLength ( 0.125 ) ;
      }   
  }

  void Execute ( const itk::Object *caller, const itk::EventObject & event ) override
  {
    // std::cout << "HelloWorld" << std::endl ;
  }

};

int main (int argc, char **argv)
{
  if ( argc <= 3 )
    {
      std::cout << "Usage: " << argv[0] << " " << "inputMovingFileName inputFixedFileName inputlabelFileName outputIaff outputLaff" << std::endl ;
      return 0 ;
    }
  
  // declare our image reader type
  typedef itk::ImageFileReader < myImageType > myFileReaderType ;

  // read the file in
  myFileReaderType::Pointer movingReader = myFileReaderType::New() ;
  movingReader->SetFileName ( argv[1] ) ;
  movingReader->Update() ;

  myFileReaderType::Pointer fixedReader = myFileReaderType::New() ;
  fixedReader->SetFileName ( argv[2] ) ;
  fixedReader->Update() ;

  myFileReaderType::Pointer labelReader = myFileReaderType::New() ;
  labelReader->SetFileName ( argv[3] ) ;
  labelReader->Update() ;

  // the registration algorithm

  // all the typedefs
  //typedef itk::ImageRegistrationMethod < myImageType, myImageType > RegistrationMethod ;
  typedef itk::MultiResolutionImageRegistrationMethod < myImageType, myImageType > RegistrationMethod ;
  typedef itk::AffineTransform < double, 3 > TransformType ;
  typedef itk::LinearInterpolateImageFunction < myImageType, double > InterpolatorType ;
  typedef itk::MeanSquaresImageToImageMetric < myImageType, myImageType > MetricType ;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType ;
  
  // all the pointers
  RegistrationMethod::Pointer regMethod = RegistrationMethod::New() ;
  TransformType::Pointer transform = TransformType::New() ;
  InterpolatorType::Pointer interpolator = InterpolatorType::New() ;
  MetricType::Pointer metric = MetricType::New() ;
  OptimizerType::Pointer optimizer = OptimizerType::New() ;
  
  // parameter setup
  regMethod->SetMovingImage ( movingReader->GetOutput() ) ;
  regMethod->SetFixedImage ( fixedReader->GetOutput() ) ;
  regMethod->SetTransform ( transform ) ;
  regMethod->SetInterpolator ( interpolator ) ;
  regMethod->SetMetric ( metric ) ;
  regMethod->SetOptimizer ( optimizer ) ;

  // additions for the multi-res version
  regMethod->SetNumberOfLevels ( 3 ) ;
  regMethod->SetFixedImageRegion ( fixedReader->GetOutput()->GetLargestPossibleRegion() ) ;
  
  // look at transform params that might need to be set up
  transform->SetIdentity() ;
  regMethod->SetInitialTransformParameters ( transform->GetParameters() ) ;
  
  // look at interpolator params that might need to be set up - looks ok to leave alone

  // todo: look at metric params that might need to be set up
  //metric->SetTransform ( transform ) ;
  //metric->Initialize() ;
  // leaving alone for now - the reg method is taking care of these for us when we call update
  
  // todo: look at optimizer params that might need to be set up
  optimizer->MinimizeOn() ;
  optimizer->SetNumberOfIterations ( 40 ) ;
  std::cout << "Min: " << optimizer->GetMinimumStepLength () << std::endl ;
  std::cout << "Max: " << optimizer->GetMaximumStepLength () << std::endl ;
  optimizer->SetMinimumStepLength ( 0 ) ;
  optimizer->SetMaximumStepLength ( 0.125 ) ;

  // set up an observer
  MyObserverCommandClass::Pointer customCommand = MyObserverCommandClass::New() ;
  optimizer->AddObserver ( itk::IterationEvent (), customCommand ) ;

  // set up second observer, this time for the registration method
  MyRegistrationObserverCommandClass::Pointer regObserver = MyRegistrationObserverCommandClass::New() ;
  regMethod->AddObserver ( itk::IterationEvent(), regObserver ) ;
  
  // run the registration
  regMethod->Update() ; 
  std::cout << "GradMagTolerance: " << optimizer->GetGradientMagnitudeTolerance() << std::endl ;
  std::cout << "Stop condition: " << optimizer->GetStopCondition() << std::endl ;
  
  // get the transform from regmethod and apply it to moving image
  typedef itk::ResampleImageFilter < myImageType, myImageType > ResampleFilterType ;
  ResampleFilterType::Pointer resampleMoveFilter = ResampleFilterType::New() ;
  ResampleFilterType::Pointer resampleLabelFilter = ResampleFilterType::New() ;

  // moving 
  resampleMoveFilter->SetInput ( movingReader->GetOutput() ) ;
  resampleMoveFilter->SetReferenceImage ( fixedReader->GetOutput() ) ;
  resampleMoveFilter->UseReferenceImageOn() ;
  resampleMoveFilter->SetDefaultPixelValue ( 0 ) ;

  // label
  resampleLabelFilter->SetInput ( labelReader->GetOutput() ) ;
  resampleLabelFilter->SetReferenceImage ( fixedReader->GetOutput() ) ;
  resampleLabelFilter->UseReferenceImageOn() ;
  resampleLabelFilter->SetDefaultPixelValue ( 0 ) ;


  transform->SetParameters ( regMethod->GetLastTransformParameters () ) ;
  resampleMoveFilter->SetTransform ( transform ) ;
  resampleMoveFilter->Update () ;


  resampleLabelFilter->SetTransform ( transform ) ;
  resampleLabelFilter->Update () ;

  



  // Part 3 Deformably register 
  const unsigned int     SpaceDimension = 3;
  constexpr unsigned int SplineOrder = 3;

  using CoordinateRepType = double;

  using DeformTransformType =
    itk::BSplineTransform<CoordinateRepType, SpaceDimension, SplineOrder>;
   auto DeformTtransform = DeformTransformType::New();
  //auto DeformTtransform = TransformType::New();

  unsigned int numberOfGridNodesInOneDimension = 5;

  
    // Software Guide : BeginCodeSnippet
  DeformTransformType::PhysicalDimensionsType fixedPhysicalDimensions;
  DeformTransformType::MeshSizeType           meshSize;
  DeformTransformType::OriginType             fixedOrigin;

  for (unsigned int i = 0; i < SpaceDimension; ++i)
  {
    fixedOrigin[i] = fixedReader->GetOutput()->GetOrigin()[i];
    fixedPhysicalDimensions[i] =
      fixedReader->GetOutput()->GetSpacing()[i] *
      static_cast<double>(
        fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize()[i] - 1);
  }
  meshSize.Fill(numberOfGridNodesInOneDimension - SplineOrder);
 
  DeformTtransform->SetTransformDomainOrigin(fixedOrigin);
  DeformTtransform->SetTransformDomainPhysicalDimensions(fixedPhysicalDimensions);
  DeformTtransform->SetTransformDomainMeshSize(meshSize);
  DeformTtransform->SetTransformDomainDirection(fixedReader->GetOutput()->GetDirection());

  using ParametersType = DeformTransformType::ParametersType;
 
  unsigned int numberOfBSplineParameters =
    DeformTtransform->GetNumberOfParameters();
 
  using RigidTransformType = itk::VersorRigid3DTransform<double>;
  using OptimizerScalesType = OptimizerType::ScalesType;
  auto rigidTransform = RigidTransformType::New();

  using AffineTransformType = itk::AffineTransform<double, SpaceDimension>;
  auto affineTransform = AffineTransformType::New();
  OptimizerScalesType optimizerScales(
    rigidTransform->GetNumberOfParameters());
  optimizerScales =
    OptimizerScalesType(affineTransform->GetNumberOfParameters());
  optimizerScales = OptimizerScalesType(numberOfBSplineParameters);
  optimizerScales.Fill(1.0);
 
  optimizer->SetScales(optimizerScales);
 


  //part 4
  itk::TimeProbesCollectorBase   chronometer;
  itk::MemoryProbesCollectorBase memorymeter;
  OptimizerType::ParametersType finalParameters = DeformTtransform->GetParameters();
 
  // Report the time and memory taken by the registration
  chronometer.Report(std::cout);
  memorymeter.Report(std::cout);
 


  constexpr unsigned int ImageDimension = 3;
  using PixelType = short;
  using FixedImageType = itk::Image<PixelType, ImageDimension>;
  using MovingImageType = itk::Image<PixelType, ImageDimension>;
  using ResampleFilterType1 =
    itk::ResampleImageFilter<MovingImageType, FixedImageType>;
 
  auto resample = ResampleFilterType::New();
 
  // resample->SetTransform(DeformTtransform);
  // resample->SetInput ( fixedReader->GetOutput() ) ;
 
  // resample->SetSize(fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize());
  // resample->SetOutputOrigin(fixedReader->GetOutput()->GetOrigin());
  // resample->SetOutputSpacing(fixedReader->GetOutput()->GetSpacing());
  // resample->SetOutputDirection(fixedReader->GetOutput()->GetDirection());

  // resample->SetDefaultPixelValue(0);
  // resample->Update();


  resample->SetInput ( resampleLabelFilter->GetOutput() ) ;
  resample->SetReferenceImage ( fixedReader->GetOutput() ) ;
  resample->UseReferenceImageOn() ;
  resample->SetDefaultPixelValue ( 0 ) ;
  resample->SetTransform ( DeformTtransform ) ;
  resample->Update () ;


  // transform->SetParameters ( regMethod->GetLastTransformParameters () ) ;
  // resampleMoveFilter->SetTransform ( transform ) ;
  // resampleMoveFilter->Update () ;


  // resampleLabelFilter->SetTransform ( transform ) ;
  // resampleLabelFilter->Update () ;


   // part5
  using StatisticsFilterType = itk::StatisticsImageFilter<myImageType>;
  auto statistics1 = StatisticsFilterType::New();
  statistics1->SetInput(resample->GetOutput());
  statistics1-> Update();
  std::cout<<"print out"<< statistics1->GetSum() << std::endl;

  
  // save the image to a different file
  typedef itk::ImageFileWriter < myImageType > myFileWriterType ;
  myFileWriterType::Pointer myMoveFileWriter = myFileWriterType::New() ;
  myMoveFileWriter->SetFileName ( argv[4] ) ;
  myMoveFileWriter->SetInput ( resampleMoveFilter->GetOutput() ) ; 
  myMoveFileWriter->Write();

  myFileWriterType::Pointer myLabelFileWriter = myFileWriterType::New() ;
  myLabelFileWriter->SetFileName ( argv[5] ) ;
  myLabelFileWriter->SetInput ( resampleLabelFilter->GetOutput() ) ; 
  myLabelFileWriter->Write();

  myFileWriterType::Pointer myLMWriter = myFileWriterType::New() ;
  myLMWriter->SetFileName ( argv[6] ) ;
  myLMWriter->SetInput ( resample->GetOutput() ) ; 
  myLMWriter->Write();
  
  return 0 ;
}

