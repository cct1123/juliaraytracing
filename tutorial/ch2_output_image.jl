function progress_bar(i, n)
    percent = floor(Int, (i / n) * 100)  # Calculate percentage
    bar_length = 100  # Length of the bar (in characters)
    progress = floor(Int, (i / n) * bar_length)  # Calculate progress
    # Create the bar string
    bar = "|" * repeat("█", progress) * repeat(" ", bar_length - progress) * "|"
    # Print the progress bar
    print("\r", "Rendering ", bar, " ", percent, "%")
    flush(stdout)  # Ensure the output is immediately printed
    if i == n
        println("") # move to the next line
    end
end

function write_ppm(filename, width, height)
    # Open the file for writing
    open(filename, "w") do file
        # PPM header
        println(file, "P3")
        println(file, "$width $height")
        println(file, "255")
        # Loop through each pixel with progress tracking
        for j in 0:height-1
            # Print progress on the same line 
            progress_bar(j, height-1)
            for i in 0:width-1
                r = i / (width - 1)  # Red component
                g = j / (height - 1) # Green component
                b = 0.25             # Blue component (constant)

                # Scale to [0, 255]
                ir = Int(round(255.999 * r))
                ig = Int(round(255.999 * g))
                ib = Int(round(255.999 * b))

                # Write the pixel color to the file
                println(file, "$ir $ig $ib")
            end
        end
        println("\rDone.                 ")
    end
end

# Call the function to create an image
width = 256
height = 256
write_ppm("images/output.ppm", width, height)
