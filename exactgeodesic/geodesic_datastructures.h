class geodesicPoints {
public:
    /// Contructor. Takes the shapes polyData as an argument.
    geodesicPoints(const double* points, int nVertices) : n(nVertices)
      {
        points_ = points;
      }

    /// returns the size of the 1D array containing all the vertex coordinates
    int size() {
        return (3 * n);
    }

    /// returns the (i % 3)-th component of the (i / 3)-th vertex.
    double operator[](int i) {
        return points_[i];
    }

    const double* points_;
    const int n;
};

/// \brief Nested class for the efficient transfer of the face data from the
/// vtk data structures to the data structures of the geodesic library.
///
class geodesicFaces {
public:
    /// Contructor. Takes the shapes polyData as an argument.
    geodesicFaces(const double* faces, const int nFaces) : n(nFaces)
      {
        faces_ = faces;
      }

    /// Returns the size of the 1D array containing all the vertex coordinates
    int size() {
        return (3 * n);
    }

    /// Returns the (i % 3)-th vertex ID of the (i / 3)-th face.
    int operator[](int i) {
        return int(faces_[i]);
    }

    const double* faces_;
    const int n;
};
