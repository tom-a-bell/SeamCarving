/**
 * @author Tom
 *
 */
public class SeamCarver
{
    /*
     * Internal deep copy of the picture to be analysed.
     */
    private Picture picture;

    /*
     * Image width and height.
     */
    private int width, height;

    /*
     * Pixel class.
     */
    private static final class Pixel implements Comparable<Pixel>
    {
        /**
         * The x- and y-coordinates of the pixel.
         */
        private final int x, y;

        /**
         * Creates the pixel (x, y).
         * @param x coordinate
         * @param y coordinate
         */
        public Pixel(final int x, final int y)
        {
            this.x = x;
            this.y = y;
        }

        /**
         * Determines if this pixel is lexicographically smaller than that pixel
         * by comparing y-coordinates and breaking ties by x-coordinates.
         * @param that pixel
         * @return -1 if less than, 0 if equal to, +1 if greater than that point.
         */
        public int compareTo(final Pixel that)
        {
            if (this.y < that.y) return -1;
            if (this.y > that.y) return +1;
            if (this.x < that.x) return -1;
            if (this.x > that.x) return +1;
            return 0;
        }

        /**
         * Returns a string representation of the pixel.
         * @return the pixel coordinates as a string.
         */
        public String toString()
        {
            return "(" + x + "," + y + ")";
        }
    }

    /*
     * Array of energy values for each pixel.
     */
    private double[][] energy;

    /*
     * Keep track of which pixels have been examined in the depth first search.
     */
    private boolean[][] marked;

    /*
     * Topological order of the pixels returned by the depth first search.
     */
    private Stack<Pixel> topologicalOrder;

    /*
     * Stored locations of the "phantom" top, bottom, left and right pixels.
     */
    private Pixel top, bottom, left, right;
    
    /*
     * Distance to each pixel returned by the shortest path algorithm.
     */
    private double[][] distTo;

    /*
     * Connected pixels in the horizontal or vertical seam.
     */
    private ST<Pixel, Pixel> pixelTo;

    /**
     * Default constructor, creating a deep copy of the supplied image.
     * @param picture - the picture to analyse
     */
    public SeamCarver(final Picture picture)
    {
        this.picture = new Picture(picture);
        this.width   = picture.width();
        this.height  = picture.height();

        // Define the positions of the phantom pixels connected to
        // the top and bottom rows and the left and right columns
        top = new Pixel(width, 0);
        bottom = new Pixel(width, height);
        left = new Pixel(0, height);
        right = new Pixel(width, height);
    }

    /**
     * Returns the current picture.
     * @return picture
     */
    public Picture picture()
    { return picture; }

    /**
     * Returns the width of the current picture.
     * @return width
     */
    public int width()
    { return width; }

    /**
     * Returns the height of the current picture.
     * @return height
     */
    public int height()
    { return height; }

    /**
     * Calculates the energy of the pixel at column x and row y.
     * @param x pixel column
     * @param y pixel row
     * @return energy
     */
    public double energy(final int x, final int y)
    {
        if (x == width || y == height)
             return 0;
        if ((x > 0) && (x < width - 1) && (y > 0) && (y < height - 1))
             return calculateEnergy(x + 1, y, x - 1, y)
                  + calculateEnergy(x, y + 1, x, y - 1);
        else return 195075;
    }

    /**
     * Returns the sequence of indices for a horizontal seam.
     * @return array of pixel indices in the seam
     */
    public int[] findHorizontalSeam()
    {
        if (energy == null)
        {
            energy = new double[width + 1][height + 1];
            for (int x = 0; x <= width; x++)
                for (int y = 0; y <= height; y++)
                    energy[x][y] = energy(x, y);
        }

        // Obtain the topological sort order by performing a depth first search,
        // starting with the left "phantom pixel" that connects to all pixels in
        // the left column and ending with the right "phantom pixel" that connects
        // to all pixels in the right column.
        marked = new boolean[width + 1][height + 1];
        topologicalOrder = new Stack<Pixel>();

        depthFirstSearch(left, true);
        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++)
                if (!marked[x][y]) depthFirstSearch(new Pixel(x, y), true);

        // Find the shortest vertical path using the topological sort order
        distTo = new double[width + 1][height + 1];
        pixelTo = new ST<Pixel, Pixel>();

        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++)
                distTo[x][y] = Double.POSITIVE_INFINITY;
        distTo[right.x][right.y] = Double.POSITIVE_INFINITY;
        distTo[left.x][left.y] = 0;

        for (Pixel p : topologicalOrder)
            for (Pixel q : horizontalNeighbours(p))
                relax(p, q);

        // Return the array of row indices in the shortest horizontal seam by
        // tracing back the pixelTo connections from the right "phantom pixel"
        int[]seam = new int[width];
        for (Pixel p = pixelTo.get(right); pixelTo.contains(p); p = pixelTo.get(p))
            seam[p.x] = p.y;

        return seam;
    }

    /**
     * Returns the sequence of indices for a vertical seam.
     * @return array of pixel indices in the seam
     */
    public int[] findVerticalSeam()
    {
        if (energy == null)
        {
            energy = new double[width + 1][height + 1];
            for (int x = 0; x <= width; x++)
                for (int y = 0; y <= height; y++)
                    energy[x][y] = energy(x, y);
        }

        // Obtain the topological sort order by performing a depth first search,
        // starting with the top "phantom pixel" that connects to all pixels in
        // the top row and ending with the bottom "phantom pixel" that connects
        // to all pixels in the bottom row.
        marked = new boolean[width + 1][height + 1];
        topologicalOrder = new Stack<Pixel>();

        depthFirstSearch(top, false);
        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++)
                if (!marked[x][y]) depthFirstSearch(new Pixel(x, y), false);

        // Find the shortest vertical path using the topological sort order
        distTo = new double[width + 1][height + 1];
        pixelTo = new ST<Pixel, Pixel>();

        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++)
                distTo[x][y] = Double.POSITIVE_INFINITY;
        distTo[bottom.x][bottom.y] = Double.POSITIVE_INFINITY;
        distTo[top.x][top.y] = 0;

        for (Pixel p : topologicalOrder)
            for (Pixel q : verticalNeighbours(p))
                relax(p, q);

        // Return the array of column indices in the shortest vertical seam by
        // tracing back the pixelTo connections from the bottom "phantom pixel"
        int[]seam = new int[height];
        for (Pixel p = pixelTo.get(bottom); pixelTo.contains(p); p = pixelTo.get(p))
            seam[p.y] = p.x;

        return seam;
    }

    /**
     * Removes the horizontal seam from the picture.
     * @param a array of pixel indices in the seam
     */
    public void removeHorizontalSeam(final int[] a)
    {
    }

    /**
     * Removes the vertical seam from the picture.
     * @param a array of pixel indices in the seam
     */
    public void removeVerticalSeam(final int[] a)
    {
    }

    private double calculateEnergy(int x1, int y1, int x2, int y2)
    {
        int R = picture.get(x1, y1).getRed()   - picture.get(x2, y2).getRed();
        int G = picture.get(x1, y1).getGreen() - picture.get(x2, y2).getGreen();
        int B = picture.get(x1, y1).getBlue()  - picture.get(x2, y2).getBlue();
        return Math.pow(R, 2) + Math.pow(G, 2) + Math.pow(B, 2);
    }

    /*
     * Performs a recursive depth first search of all connected pixels in the 
     * specified direction (horizontal = <True>, vertical = <False>).
     */
    private void depthFirstSearch(Pixel p, boolean horizontal)
    {
        marked[p.x][p.y] = true;

        // Find the connected pixels (either horizontally or vertically
        Stack<Pixel> connectedPixels = new Stack<Pixel>();
        if (horizontal) connectedPixels = horizontalNeighbours(p);
        else            connectedPixels = verticalNeighbours(p);

        // Traverse each connected pixel
        for (Pixel q : connectedPixels)
            if (!marked[q.x][q.y]) depthFirstSearch(q, horizontal);

        topologicalOrder.push(p);
    }

    private void relax(Pixel p, Pixel q)
    {
        if (distTo[q.x][q.y] > distTo[p.x][p.y] + energy[q.x][q.y])
        {
            distTo[q.x][q.y] = distTo[p.x][p.y] + energy[q.x][q.y];
            pixelTo.put(q, p);
        }
    }

    private Stack<Pixel> horizontalNeighbours(Pixel p)
    {
        Stack<Pixel> horizontalNeighbours = new Stack<Pixel>();
        if (p == left)
            for (int y = 0; y < height; y++)
                horizontalNeighbours.push(new Pixel(0, y));
        else if (p.x == width - 1)
            horizontalNeighbours.push(right);
        else if (p.x < width - 1)
        {
            if (p.y > 0)
                horizontalNeighbours.push(new Pixel(p.x + 1, p.y - 1));
            horizontalNeighbours.push(new Pixel(p.x + 1, p.y));
            if (p.y < height - 1)
                horizontalNeighbours.push(new Pixel(p.x + 1, p.y + 1));
        }
        return horizontalNeighbours;
    }

    private Stack<Pixel> verticalNeighbours(Pixel p)
    {
        Stack<Pixel> verticalNeighbours = new Stack<Pixel>();
        if (p == top)
            for (int x = 0; x < width; x++)
                verticalNeighbours.push(new Pixel(x, 0));
        else if (p.y == height - 1)
            verticalNeighbours.push(bottom);
        else if (p.y < height - 1)
        {
            if (p.x > 0)
                verticalNeighbours.push(new Pixel(p.x - 1, p.y + 1));
            verticalNeighbours.push(new Pixel(p.x, p.y + 1));
            if (p.x < width - 1)
                verticalNeighbours.push(new Pixel(p.x + 1, p.y + 1));
        }
        return verticalNeighbours;
    }

    /**
     * Unit tests.
     * @param args command line arguments
     */
    public static void main(final String[] args)
    {
    }

}
