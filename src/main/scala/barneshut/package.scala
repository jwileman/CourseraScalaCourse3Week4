import common._
import barneshut.conctrees._

package object barneshut {

  class Boundaries {
    var minX = Float.MaxValue

    var minY = Float.MaxValue

    var maxX = Float.MinValue

    var maxY = Float.MinValue

    def width = maxX - minX

    def height = maxY - minY

    def size = math.max(width, height)

    def centerX = minX + width / 2

    def centerY = minY + height / 2

    override def toString = s"Boundaries($minX, $minY, $maxX, $maxY)"
  }

  sealed abstract class Quad {
    def massX: Float

    def massY: Float

    def mass: Float

    def centerX: Float

    def centerY: Float

    def size: Float

    def total: Int

    def insert(b: Body): Quad
  }

  case class Empty(centerX: Float, centerY: Float, size: Float) extends Quad {
    def massX: Float = centerX
    def massY: Float = centerY
    def mass: Float = 0
    def total: Int = 0
    def insert(b: Body): Quad = Leaf(centerX, centerY, size, Seq(b))
  }

  case class Fork(
    nw: Quad, ne: Quad, sw: Quad, se: Quad
  ) extends Quad {
    val centerX: Float = nw.centerX + nw.size / 2
    val centerY: Float = nw.centerY + nw.size / 2
    val size: Float = nw.size * 2
    val mass: Float = nw.mass + ne.mass + sw.mass + se.mass
    val massX: Float = if (mass == 0) centerX else (nw.mass * nw.massX + ne.mass * ne.massX + sw.mass * sw.massX + se.mass * se.massX) / mass
    val massY: Float = if (mass == 0) centerY else (nw.mass * nw.massY + ne.mass * ne.massY + sw.mass * sw.massY + se.mass * se.massY) / mass
    val total: Int = nw.total + ne.total + sw.total + se.total

    def insert(b: Body): Fork = {
      val halfSize = size / 2
      val left = centerX - halfSize <= b.x && b.x < centerX
      val right = centerX <= b.x && b.x <= centerX + halfSize
      val top = centerY - halfSize <= b.y && b.y < centerY
      val bottom = centerY <= b.y && b.y <= centerY + halfSize

      if (left && top) {
        //nw
        Fork(nw.insert(b), ne, sw, se)
      } else if (right && top) {
        //ne
        Fork(nw, ne.insert(b), sw, se)
      } else if (left && bottom) {
        //sw
        Fork(nw, ne, sw.insert(b), se)
      } else if (right && bottom) {
        //se
        Fork(nw, ne, sw, se.insert(b))
      } else {
        this
      }
    }
  }

  case class Leaf(centerX: Float, centerY: Float, size: Float, bodies: Seq[Body])
  extends Quad {
    val (mass, massX, massY) = (bodies.foldLeft(0f)(_ + _.mass) : Float, bodies.foldLeft(0f)((sum, b) => sum + b.mass * b.x) / bodies.foldLeft(0f)(_ + _.mass) : Float, bodies.foldLeft(0f)((sum, b) => sum + b.mass * b.y) / bodies.foldLeft(0f)(_ + _.mass)  : Float)
    val total: Int = bodies.size
    def insert(b: Body): Quad = {
      if (size > minimumSize) {
        val quarterSize = size / 4
        val halfSize = size / 2
        val emptyQuadtree = Fork(Empty(centerX - quarterSize, centerY - quarterSize, halfSize), Empty(centerX + quarterSize, centerY - quarterSize, halfSize), Empty(centerX - quarterSize, centerY + quarterSize, halfSize), Empty(centerX + quarterSize, centerY + quarterSize, halfSize))
        bodies.foldLeft(emptyQuadtree)(_.insert(_))
      } else {
        Leaf(centerX, centerY, size, bodies :+ b)
      }
    }
  }

  def minimumSize = 0.00001f

  def gee: Float = 100.0f

  def delta: Float = 0.01f

  def theta = 0.5f

  def eliminationThreshold = 0.5f

  def force(m1: Float, m2: Float, dist: Float): Float = gee * m1 * m2 / (dist * dist)

  def distance(x0: Float, y0: Float, x1: Float, y1: Float): Float = {
    math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)).toFloat
  }

  class Body(val mass: Float, val x: Float, val y: Float, val xspeed: Float, val yspeed: Float) {

    def updated(quad: Quad): Body = {
      var netforcex = 0.0f
      var netforcey = 0.0f

      def addForce(thatMass: Float, thatMassX: Float, thatMassY: Float): Unit = {
        val dist = distance(thatMassX, thatMassY, x, y)
        /* If the distance is smaller than 1f, we enter the realm of close
         * body interactions. Since we do not model them in this simplistic
         * implementation, bodies at extreme proximities get a huge acceleration,
         * and are catapulted from each other's gravitational pull at extreme
         * velocities (something like this:
         * http://en.wikipedia.org/wiki/Interplanetary_spaceflight#Gravitational_slingshot).
         * To decrease the effect of this gravitational slingshot, as a very
         * simple approximation, we ignore gravity at extreme proximities.
         */
        if (dist > 1f) {
          val dforce = force(mass, thatMass, dist)
          val xn = (thatMassX - x) / dist
          val yn = (thatMassY - y) / dist
          val dforcex = dforce * xn
          val dforcey = dforce * yn
          netforcex += dforcex
          netforcey += dforcey
        }
      }

      def traverse(quad: Quad): Unit = (quad: Quad) match {
        case Empty(_, _, _) =>
          // no force
        case Leaf(_, _, _, bodies) =>
          // add force contribution of each body by calling addForce
          bodies.foreach(body => addForce(body.mass, body.x, body.y))
        case Fork(nw, ne, sw, se) =>
          // see if node is far enough from the body,
          // or recursion is needed

          // distance between the center of mass of quad and the particle
          val dist = distance(x, y, quad.massX, quad.massY)

          if (quad.size / dist < theta) {
            // a fork quadtree that is sufficiently far away acts as a single point of mass
            addForce(quad.mass, quad.massX, quad.massY)
          } else {
            // a fork quadtree that is not sufficiently far away must be recursively traversed
            traverse(nw)
            traverse(ne)
            traverse(sw)
            traverse(se)
          }
      }

      traverse(quad)

      val nx = x + xspeed * delta
      val ny = y + yspeed * delta
      val nxspeed = xspeed + netforcex / mass * delta
      val nyspeed = yspeed + netforcey / mass * delta

      new Body(mass, nx, ny, nxspeed, nyspeed)
    }

  }

  val SECTOR_PRECISION = 8

  class SectorMatrix(val boundaries: Boundaries, val

  sectorPrecision: Int) {
    val sectorSize = boundaries.size / sectorPrecision
    val matrix = new Array[ConcBuffer[Body]](sectorPrecision * sectorPrecision)
    for (i <- 0 until matrix.length) matrix(i) = new ConcBuffer

    def +=(b: Body): SectorMatrix = {
      val leftX = boundaries.centerX - boundaries.size / 2
      val rightX = boundaries.centerX + boundaries.size / 2
      val topY = boundaries.centerY - boundaries.size / 2
      val bottomY = boundaries.centerY + boundaries.size / 2

      val correctX = if (b.x < leftX) leftX else if (b.x > rightX) rightX else b.x
      val col = Math.min(((correctX - leftX) / sectorSize).toInt, sectorPrecision - 1)

      val correctY = if (b.y < topY) topY else if (b.y > bottomY) bottomY else b.y
      val row = Math.min(((correctY - topY) / sectorSize).toInt, sectorPrecision - 1)

      this(col, row) += b

      this
    }

    def apply(col: Int, row: Int) = matrix(row * sectorPrecision + col)

    def combine(that: SectorMatrix): SectorMatrix = {
      //for (i <- this.matrix.size) yield this.matrix(i) c
      val newSectorMatrix = new SectorMatrix(this.boundaries, this.sectorPrecision)

      for (col <- 0 until sectorPrecision; row <- 0 until sectorPrecision) {
        this(col, row).combine(that(col, row)).foreach(newSectorMatrix += _)
      }

      newSectorMatrix
    }

    /**
      * JW note to self: achievedParallelism is 1 (thread), 4 (threads), 16 (threads), 64 (threads) etc.
      */
    def toQuad(parallelism: Int): Quad = {
      def BALANCING_FACTOR = 4
      def quad(col: Int, row: Int, span: Int, achievedParallelism: Int): Quad = {
        if (span == 1) {
          val sectorSize = boundaries.size / sectorPrecision
          val centerX = boundaries.minX + col * sectorSize + sectorSize / 2
          val centerY = boundaries.minY + row * sectorSize + sectorSize / 2
          var emptyQuad: Quad = Empty(centerX, centerY, sectorSize)
          val sectorBodies = this(col, row)
          sectorBodies.foldLeft(emptyQuad)(_ insert _)
        } else {
          val nspan = span / 2
          val nAchievedParallelism = achievedParallelism * 4
          val (nw, ne, sw, se) =
            if (parallelism > 1 && achievedParallelism < parallelism * BALANCING_FACTOR) parallel(
              quad(col, row, nspan, nAchievedParallelism),
              quad(col + nspan, row, nspan, nAchievedParallelism),
              quad(col, row + nspan, nspan, nAchievedParallelism),
              quad(col + nspan, row + nspan, nspan, nAchievedParallelism)
            ) else (
              quad(col, row, nspan, nAchievedParallelism),
              quad(col + nspan, row, nspan, nAchievedParallelism),
              quad(col, row + nspan, nspan, nAchievedParallelism),
              quad(col + nspan, row + nspan, nspan, nAchievedParallelism)
            )
          Fork(nw, ne, sw, se)
        }
      }

      quad(0, 0, sectorPrecision, 1)
    }

    override def toString = s"SectorMatrix(#bodies: ${matrix.map(_.size).sum})"
  }

  class TimeStatistics {
    private val timeMap = collection.mutable.Map[String, (Double, Int)]()

    def clear() = timeMap.clear()

    def timed[T](title: String)(body: =>T): T = {
      var res: T = null.asInstanceOf[T]
      val totalTime = /*measure*/ {
        val startTime = System.currentTimeMillis()
        res = body
        (System.currentTimeMillis() - startTime)
      }

      timeMap.get(title) match {
        case Some((total, num)) => timeMap(title) = (total + totalTime, num + 1)
        case None => timeMap(title) = (0.0, 0)
      }

      println(s"$title: ${totalTime} ms; avg: ${timeMap(title)._1 / timeMap(title)._2}")
      res
    }

    override def toString = {
      timeMap map {
        case (k, (total, num)) => k + ": " + (total / num * 100).toInt / 100.0 + " ms"
      } mkString("\n")
    }
  }
}
