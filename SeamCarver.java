/* *****************************************************************************
 *  Name:
 *  Date:
 *  Description:
 **************************************************************************** */

import edu.princeton.cs.algs4.Picture;

public class SeamCarver {

    private Picture picture;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture)
    {
        if (picture == null) {
            throw new IllegalArgumentException();
        }
        // The data type may not mutate the Picture argument to the constructor
        this.picture = new Picture(picture);
    }

    // current picture
    public Picture picture()
    {
        return new Picture(picture);
    }

    // width of current picture
    public int width()
    {
        return picture.width();
    }

    // height of current picture
    public int height()
    {
        return picture.height();
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y)
    {
        if (x < 0 || x > width() - 1 || y < 0 || y > height() - 1) {
            throw new IllegalArgumentException();
        }

        if (x == 0 || x == width() - 1 || y == 0 || y == height() - 1) {
            return 1000;
        }

        int rgb_x = picture.getRGB(x - 1, y);
        int r_x = (rgb_x >> 16) & 0xFF, g_x = (rgb_x >> 8) & 0xFF, b_x = (rgb_x >> 0) & 0xFF;
        int rgb_x1 = picture.getRGB(x + 1, y);
        int r_x1 = (rgb_x1 >> 16) & 0xFF, g_x1 = (rgb_x1 >> 8) & 0xFF, b_x1 = (rgb_x1 >> 0) & 0xFF;
        int rgb_y = picture.getRGB(x, y - 1);
        int r_y = (rgb_y >> 16) & 0xFF, g_y = (rgb_y >> 8) & 0xFF, b_y = (rgb_y >> 0) & 0xFF;
        int rgb_y1 = picture.getRGB(x, y + 1);
        int r_y1 = (rgb_y1 >> 16) & 0xFF, g_y1 = (rgb_y1 >> 8) & 0xFF, b_y1 = (rgb_y1 >> 0) & 0xFF;

        double grad_Rx = r_x - r_x1;
        double grad_Gx = g_x - g_x1;
        double grad_Bx = b_x - b_x1;

        double grad_Ry = r_y - r_y1;
        double grad_Gy = g_y - g_y1;
        double grad_By = b_y - b_y1;

        return Math.sqrt(grad_Rx * grad_Rx + grad_Gx * grad_Gx + grad_Bx * grad_Bx +
                grad_Ry * grad_Ry + grad_Gy * grad_Gy + grad_By * grad_By);
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam()
    {
        transposePicture();
        int[] seam = findVerticalSeam();
        transposePicture();

        return seam;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam()
    {
        // cache the energy levels for each pixel
        double[][] energy = new double[picture.height()][];
        // the lowest energy cost to each pixel
        double[][] distTo = new double[picture.height()][];
        // the lowest energy path to each pixel
        int[][] edgeTo = new int[picture.height()][];

        for (int row = 0; row < picture.height(); row++) {
            energy[row] = new double[picture.width()];
            distTo[row] = new double[picture.width()];
            edgeTo[row] = new int[picture.width()];
            for (int col = 0; col < picture.width(); col++) {
                energy[row][col] = energy(col, row);
                // the distance to any of the pixels in the first row is zero
                distTo[row][col] = row == 0 ? 0 : Double.POSITIVE_INFINITY;
                edgeTo[row][col] = -1;
            }
        }

        // traverse image in topological order (skip last row, it has no children)
        for (int row = 0; row < picture.height() - 1; row++) {
            for (int col = 0; col < picture.width(); col++) {
                // relax each edge
                int next_row = row + 1;
                for (int next_col = col - 1; next_col <= col + 1; next_col++) {
                    if (next_col >= 0 && next_col < picture.width()) {
                        double weight = energy[next_row][next_col];
                        if (distTo[next_row][next_col] > distTo[row][col] + weight) {
                            distTo[next_row][next_col] = distTo[row][col] + weight;
                            // save the column in the edge-to graph (row is implicit)
                            edgeTo[next_row][next_col] = col;
                        }
                    }
                }
            }
        }

        // the lowest energy vertical path through the image
        int[] seam = new int[picture.height()];

        // find the lowest-energy path that reaches the last row
        double min = Double.POSITIVE_INFINITY;
        int low_col = -1;
        for (int col = 0; col < picture.width(); col++) {
            if (min > distTo[picture.height() - 1][col]) {
                low_col = col;
                min = distTo[picture.height() - 1][col];
            }
        }
        assert(low_col > 0);
        for (int row = picture.height() - 1; row >= 0; row--) {
            seam[row] = low_col;
            low_col = edgeTo[row][low_col];
        }

        return seam;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam)
    {
        transposePicture();
        removeVerticalSeam(seam);
        transposePicture();
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam)
    {
        if (picture.width() <= 1) {
            throw new IllegalArgumentException("insufficient picture width");
        }
        validateSeam(seam);

        // replace picture with smaller image - very slow
        Picture trimmed = new Picture(picture.width() - 1, picture.height());
        for (int row = 0; row < picture.height(); row++) {
            int trimColumn = seam[row];
            for (int col = 0; col < trimColumn; col++) {
                trimmed.setRGB(col, row, picture.getRGB(col, row));
            }
            for (int col = trimColumn + 1; col < picture.width(); col++) {
                trimmed.setRGB(col - 1, row, picture.getRGB(col, row));
            }
        }
        this.picture = trimmed;
    }

    private void transposePicture()
    {
        // swap width and height for new image
        Picture transpose = new Picture(picture.height(), picture.width());
        for (int col = 0; col < picture.width(); col++) {
            for (int row = 0; row < picture.height(); row++) {
                transpose.setRGB(row, col, picture.getRGB(col, row));
            }
        }
        this.picture = transpose;
    }

    private void validateSeam(int[] seam)
    {
        if (seam == null) {
            throw new IllegalArgumentException("seam is null");
        }

        if (seam.length != picture.height()) {
            throw new IllegalArgumentException("seam length does not match range of " + picture.height());
        }

        int prev = seam[0];
        if (prev < 0 || prev >= picture.width()) {
            throw new IllegalArgumentException("seam entry " + prev + " is outside prescribed range of " + picture.width());
        }
        for (int i = 1; i < seam.length; ++i) {
            int curr = seam[i];
            if (Math.abs(curr - prev) > 1) {
                throw new IllegalArgumentException("adjacent seam entries are too var apart");
            }
            if (curr < 0 || curr >= picture.width()) {
                throw new IllegalArgumentException("seam entry " + i + " is outside prescribed range of " + picture.width());
            }
            prev = curr;
        }
    }
}
