using Images, FileIO

function ppm_to_jpg(ppm_filename::String, jpg_filename::String)
    # Read the PPM file
    img = load(ppm_filename)

    # Save the image as a JPEG file
    save(jpg_filename, img)
end