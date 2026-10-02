#include <Rcpp.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNiftiImageIOFactory.h>
#include <itkTransformFactoryBase.h>
#include <itkHDF5TransformIOFactory.h>
#include <itkMatlabTransformIOFactory.h>
#include <itkTxtTransformIOFactory.h>
#include <itkTransformFileReader.h>
#include <itkCompositeTransform.h>
#include <itkIdentityTransform.h>
#include <itkResampleImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkCastImageFilter.h>
#include <cmath>

namespace {
using Image = itk::Image<double, 3>;
using Transform = itk::Transform<double, 3, 3>;
using Composite = itk::CompositeTransform<double, 3>;

void register_io() {
  static const bool registered = []() {
    itk::NiftiImageIOFactory::RegisterOneFactory();
    itk::HDF5TransformIOFactory::RegisterOneFactory();
    itk::MatlabTransformIOFactory::RegisterOneFactory();
    itk::TxtTransformIOFactory::RegisterOneFactory();
    itk::TransformFactoryBase::RegisterDefaultTransforms();
    return true;
  }();
  (void) registered;
}
Image::Pointer read_image(const std::string& path) {
  auto reader = itk::ImageFileReader<Image>::New();
  reader->SetFileName(path);
  reader->Update();
  Image::Pointer image = reader->GetOutput();
  image->DisconnectPipeline();
  return image;
}
void write_image(Image* image, const std::string& path, bool singleprecision) {
  if (singleprecision) {
    using FloatImage = itk::Image<float, 3>;
    auto cast = itk::CastImageFilter<Image, FloatImage>::New();
    cast->SetInput(image);
    auto writer = itk::ImageFileWriter<FloatImage>::New();
    writer->SetFileName(path);
    writer->SetInput(cast->GetOutput());
    writer->SetUseCompression(true);
    writer->Update();
  } else {
    auto writer = itk::ImageFileWriter<Image>::New();
    writer->SetFileName(path);
    writer->SetInput(image);
    writer->SetUseCompression(true);
    writer->Update();
  }
}
}

// Scalar 3-D images only. ITK uses physical LPS coordinates and zero-based
// indices internally; no manual NIfTI RAS/LPS conversion is performed.
// [[Rcpp::export]]
void qsi_apply_transforms_cpp(std::string fixed, std::string moving,
                              Rcpp::CharacterVector transformlist,
                              std::string output, std::string interpolator,
                              double defaultvalue, int nthread,
                              bool singleprecision) {
  if (nthread < 1) Rcpp::stop("nthread must be positive");
  if (!std::isfinite(defaultvalue)) Rcpp::stop("defaultvalue must be finite");
  register_io();
  try {
    auto reference = read_image(fixed);
    auto input = read_image(moving);
    auto composite = Composite::New();
    for (R_xlen_t i = 0; i < transformlist.size(); ++i) {
      Rcpp::checkUserInterrupt();
      auto reader = itk::TransformFileReaderTemplate<double>::New();
      reader->SetFileName(Rcpp::as<std::string>(transformlist[i]));
      reader->Update();
      const auto* transforms = reader->GetTransformList();
      if (transforms->empty()) Rcpp::stop("Transform file is empty");
      for (const auto& item : *transforms) {
        auto* transform = dynamic_cast<Transform*>(item.GetPointer());
        if (!transform) Rcpp::stop("Expected a three-dimensional double-precision ITK transform");
        composite->AddTransform(transform);
      }
    }
    if (composite->GetNumberOfTransforms() == 0)
      composite->AddTransform(itk::IdentityTransform<double, 3>::New());
    // ITK applies the last queued transform first, as does ANTs.
    using Resampler = itk::ResampleImageFilter<Image, Image, double, double>;
    auto resampler = Resampler::New();
    resampler->SetInput(input);
    resampler->SetReferenceImage(reference);
    resampler->UseReferenceImageOn();
    resampler->SetTransform(composite);
    resampler->SetDefaultPixelValue(defaultvalue);
    resampler->SetNumberOfWorkUnits(static_cast<unsigned int>(nthread));
    if (interpolator == "linear") {
      resampler->SetInterpolator(itk::LinearInterpolateImageFunction<Image, double>::New());
    } else if (interpolator == "nearestNeighbor") {
      resampler->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<Image, double>::New());
    } else {
      Rcpp::stop("Supported interpolators: linear, nearestNeighbor");
    }
    resampler->Update();
    Rcpp::checkUserInterrupt();
    write_image(resampler->GetOutput(), output, singleprecision);
  } catch (const itk::ExceptionObject& error) {
    Rcpp::stop("ITK transform/resampling error: %s", error.what());
  }
}

// [[Rcpp::export]]
void qsi_image_write_cpp(std::string input, std::string output,
                         bool singleprecision) {
  register_io();
  try {
    auto image = read_image(input);
    write_image(image, output, singleprecision);
  } catch (const itk::ExceptionObject& error) {
    Rcpp::stop("ITK image-write error: %s", error.what());
  }
}
