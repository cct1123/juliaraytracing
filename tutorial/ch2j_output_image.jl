using Images, FileIO, ProgressBars

function create_image(width, height)
    # Create an empty image array (height x width x 3 for RGB)
    img = Array{RGB{N0f8}}(undef, height, width)
    # Populate the image with colors
    for j in ProgressBar(1:height)
        for i in 1:width
            r = (i - 1) / (width - 1)  # Red component
            g = (j - 1) / (height - 1) # Green component
            b = 0.25                   # Blue component (constant)
            img[j, i] = RGB(r, g, b)   # Set the pixel
        end
    end
    return img
end

# Image dimensions
width = 256
height = 256

# Generate and save the image
img = create_image(width, height)
save("images/output.png", img)
println("Image saved as output.png")