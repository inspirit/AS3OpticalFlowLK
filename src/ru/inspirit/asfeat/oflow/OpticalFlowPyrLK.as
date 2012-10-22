package ru.inspirit.asfeat.oflow
{
    import apparat.asm.IncLocalInt;
    import apparat.asm.__asm;
    import apparat.asm.__cint;
    import apparat.math.FastMath;
    import apparat.math.IntMath;

    import flash.display.BitmapData;
    import flash.filters.ColorMatrixFilter;
    import flash.geom.Matrix;
    import flash.geom.Point;
    import flash.geom.Rectangle;

    public final class OpticalFlowPyrLK
    {
        public static const OPTFLOW_USE_INITIAL_FLOW:int = 4;
        public static const OPTFLOW_GET_MIN_EIGENVALS:int = 8;

        protected var _width:int;
        protected var _height:int;
        protected var _maxLevel:int;
        protected var _winSize:int;
        protected var _winArea:int;

        protected var _sharTempBuff0:Vector.<int>;
        protected var _sharTempBuff1:Vector.<int>;
        protected var _iWinBuf:Vector.<int>;
        protected var _derivIWinBuf:Vector.<int>;
        protected var _levelsDeriv:Vector.<Vector.<int>>;

        protected var _prevPyr:Vector.<Vector.<int>>;
        protected var _nextPyr:Vector.<Vector.<int>>;
        protected var _prevPyrBmp:Vector.<BitmapData>;
        protected var _nextPyrBmp:Vector.<BitmapData>;

        protected const _zero:Point = new Point(0,0);
        protected const _pyrMat:Matrix = new Matrix(0.5, 0, 0, 0.5);
        protected const _grayMat:ColorMatrixFilter = new ColorMatrixFilter([
            0, 0, 0, 0, 0,
            0, 0, 0, 0, 0,
            .2989, .587, .114, 0, 0,
            0, 0, 0, 0, 0
        ]);

        public function OpticalFlowPyrLK()
        {
            //
        }

        public function setLevelsDerivatives(derivs:Vector.<Vector.<int>>):void
        {
            var n:int = Math.min(derivs.length, _maxLevel);
            for(var i:int = 0; i < n; ++i)
            {
                _levelsDeriv[i] = derivs[i];
            }
        }
        
        public function setBitmapImagePyramids(prevPyr:Vector.<BitmapData>, nextPyr:Vector.<BitmapData>):void
        {
            var n:int = Math.min(prevPyr.length, _maxLevel);
            for (var i:int = 0; i < n; ++i)
            {
                _prevPyrBmp[i] = prevPyr[i];
                _nextPyrBmp[i] = nextPyr[i];
            }
        }
        
        public function setVectorImagePyramids(prevPyr:Vector.<Vector.<int>>, nextPyr:Vector.<Vector.<int>>):void
        {
            var n:int = Math.min(prevPyr.length, _maxLevel);
            for (var i:int = 0; i < n; ++i)
            {
                _prevPyr[i] = prevPyr[i];
                _nextPyr[i] = nextPyr[i];
            }
        }

        public function swapImagePyramids():void
        {
            var tmp:Vector.<Vector.<int>> = _prevPyr;
            _prevPyr = _nextPyr;
            _nextPyr = tmp;

            var tmp_bmp:Vector.<BitmapData> = _prevPyrBmp;
            _prevPyrBmp = _nextPyrBmp;
            _nextPyrBmp = tmp_bmp;
        }

        public function buildNextImagePyramid(bmp:BitmapData, doGrayScale:Boolean = false):void
        {
            var n:int = _maxLevel;
            var i:int = 0;
            var rect:Rectangle = bmp.rect;
            var b:BitmapData = _nextPyrBmp[i];
            var a:BitmapData;

            if(doGrayScale)
            {
                b.applyFilter(bmp, rect, _zero, _grayMat);
            }
            else b.copyPixels(bmp, rect, _zero);

            var db:Vector.<uint> = b.getVector(rect);
            var idb:Vector.<int> = _nextPyr[i];
            var sz:int = db.length;
            var j:int;
            for (j = 0; j < sz; ++j)
            {
                idb[j] = db[j] & 0xFF;
            }
            

            ++i;
            for(; i < n; ++i)
            {
                rect.width >>= 1;
                rect.height >>= 1;

                a = b;
                b = _nextPyrBmp[i];
                b.draw(a, _pyrMat, null, null, null, true);
                //_nextPyr[i] = b.getVector(rect);
                //
                db = b.getVector(rect);
                idb = _nextPyr[i];
                sz = db.length;
                for (j = 0; j < sz; ++j)
                {
                    idb[j] = db[j] & 0xFF;
                }
            }
        }

        public function calcPrevImageDerivatives():void
        {
            var n:int = _maxLevel;
            var i:int, lev_w:int, lev_h:int, deriv_w:int, deriv_h:int, dstep:int;
            var x:int, y:int;
            var x0:int, x1:int;
            var deriv:Vector.<int>;
            var img:Vector.<int>;
            var trow0:Vector.<int> = _sharTempBuff0;
            var trow1:Vector.<int> = _sharTempBuff1;

            for(i = 0; i < n; ++i)
            {
                lev_w = _width >> i;
                lev_h = _height >> i;
                deriv_w = lev_w;
                deriv_h = lev_h;
                dstep = __cint(deriv_w * 2);
                deriv = _levelsDeriv[i];
                img = _prevPyr[i];

                // calc SharDeriv
                for(y = 0; y < lev_h; ++y)
                {
                    var srow1:int = __cint(y*lev_w);
                    var srow0:int = __cint((y > 0 ? y-1 : 1)*lev_w);
                    var srow2:int = __cint((y < lev_h-1 ? y+1 : lev_h-2)*lev_w);
                    var drow:int = __cint(y*dstep);

                    // do vertical convolution
                    for(x = 0, x1 = 1; x < lev_w; )
                    {
                        var t0:int = __cint( ((img[srow0+x]) + (img[srow2+x]))*3 + (img[srow1+x])*10 );
                        var t1:int = __cint( (img[srow2+x]) - (img[srow0+x]) );
                        trow0[x1] = t0;
                        trow1[x1] = t1;
                        __asm(IncLocalInt(x),IncLocalInt(x1));
                    }

                    // make border
                    x = __cint(lev_w + 1);
                    trow0[0] = trow0[1]; trow0[x] = trow0[lev_w];
                    trow1[0] = trow1[1]; trow1[x] = trow1[lev_w];

                    // do horizontal convolution, interleave the results and store them
                    for(x = 0; x < lev_w; )
                    {
                        t0 = __cint( (trow0[x+2] - trow0[x]) );
                        t1 = __cint( ((trow1[x+2] + trow1[x])*3 + trow1[x+1]*10) );
                        deriv[drow] = t0;
                        __asm(IncLocalInt(drow));
                        deriv[drow] = t1;
                        __asm(IncLocalInt(drow),IncLocalInt(x));
                    }
                }
            }
        }

        public function setup(width:int, height:int, winSize:int, levels:int, initBitmaps:Boolean = true):void
        {
            _width = width;
            _height = height;
            _winSize = winSize;
            _winArea = winSize*winSize;
            _maxLevel = levels;

            var delta:int = width + 2;
            _sharTempBuff0 = new Vector.<int>(delta, true);
            _sharTempBuff1 = new Vector.<int>(delta, true);
            _iWinBuf = new Vector.<int>(_winArea, true);
            _derivIWinBuf = new Vector.<int>(_winArea*2, true);
            _levelsDeriv = new Vector.<Vector.<int>>(levels, true);
            _prevPyr = new Vector.<Vector.<int>>(levels, true);
            _nextPyr = new Vector.<Vector.<int>>(levels, true);
            _prevPyrBmp = new Vector.<BitmapData>(levels, true);
            _nextPyrBmp = new Vector.<BitmapData>(levels, true);
            for(var i:int = 0; i < levels; ++i)
            {
                var lw:int = width >> i;
                var lh:int = height >> i;
                var dw:int = lw;
                var dh:int = lh;
                var dstep:int = dw * 2;
                _levelsDeriv[i] = new Vector.<int>(dstep*dh, true);

                // TODO: use BitmapData.getPixelsToVector
                _prevPyr[i] = new Vector.<int>(lw*lh, true);
                _nextPyr[i] = new Vector.<int>(lw*lh, true);
                if(initBitmaps) {
                    _prevPyrBmp[i] = new BitmapData(lw, lh, false, 0x0);
                    _nextPyrBmp[i] = new BitmapData(lw, lh, false, 0x0);
                }
            }
        }

        public function run(prevPts:Vector.<Number>, nextPts:Vector.<Number>, count:int,
                            status:Vector.<int> = null, error:Vector.<Number> = null,
                            maxIterations:int = 30, epsilon:Number = 0.01, minEigThreshold:Number = 0.0001, flags:int = 0):void
        {
            if(count <= 0) return;

            var i:int, j:int, x:int, y:int, level:int, ptidx:int;
            var lev_w:int, lev_h:int, deriv_w:int, deriv_h:int, dstep:int;
            //var deriv_off:int;
            var deriv:Vector.<int>;
            var img_prev:Vector.<int>;
            var img_next:Vector.<int>;

            var winSize:int = _winSize;
            var halfWin:Number = (winSize-1)*0.5;
            var winArea:int = winSize * winSize;
            var winArea2:int = (2 * winArea);
            var invWin:Number = 1.0 / Number(winArea);
            var IWinBuf:Vector.<int> = _iWinBuf;
            var derivIWinBuf:Vector.<int> = _derivIWinBuf;
            var IWinBuf_step:int = winSize;
            var derivIWinBuf_step:int = winSize*2;
            var lev_sc:Number;
            var prevPt_x:Number, prevPt_y:Number, nextPt_x:Number, nextPt_y:Number;
            var nextPoint_x:Number, nextPoint_y:Number;
            var prevDelta_x:Number, prevDelta_y:Number;
            var iprevPt_x:int, iprevPt_y:int, inextPt_x:int, inextPt_y:int;
            var inextPoint_x:int, inextPoint_y:int;

            // squared epsilon is used
            const FLT_EPSILON:Number = 0.00000011920929;
            epsilon *= epsilon;
            var _status:Vector.<int>;
            if(null == status) status = _status = new Vector.<int>(count);
            for(i = 0; i < count; ++i) status[i] = 1;

            var maxLevel:int = _maxLevel - 1;
            level = maxLevel;

            for(; level >= 0; --level)
            {
                lev_sc = (1./(1 << level));
                lev_w = _width >> level;
                lev_h = _height >> level;
                deriv_w = lev_w;
                deriv_h = lev_h;
                dstep = __cint(deriv_w * 2);
                deriv = _levelsDeriv[level];
                img_prev = _prevPyr[level];
                img_next = _nextPyr[level];
                
                var brd_l:int = 0;
                var brd_t:int = 0;
                var brd_r:int = lev_w - winSize;
                var brd_b:int = lev_h - winSize;

                // iterate through points
                for(ptidx = 0; ptidx < count; ++ptidx)
                {
                    const id0:int = ptidx << 1;
                    const id1:int = __cint(id0 + 1);
                    prevPt_x = prevPts[id0]*lev_sc;
                    prevPt_y = prevPts[id1]*lev_sc;

                    if( level == maxLevel )
                    {
                        if( flags & OPTFLOW_USE_INITIAL_FLOW )
                        {
                            nextPt_x = nextPts[id0]*lev_sc;
                            nextPt_y = nextPts[id1]*lev_sc;
                        } else {
                            nextPt_x = prevPt_x;
                            nextPt_y = prevPt_y;
                        }
                    }
                    else {
                        nextPt_x = nextPts[id0]*2.0;
                        nextPt_y = nextPts[id1]*2.0;
                    }
                    nextPts[id0] = nextPt_x;
                    nextPts[id1] = nextPt_y;

                    prevPt_x -= halfWin;
                    prevPt_y -= halfWin;
                    iprevPt_x = (prevPt_x);
                    iprevPt_y = (prevPt_y);
                    
                    x = int(iprevPt_x <= brd_l) 
                        | int(iprevPt_x >= brd_r) 
                        | int(iprevPt_y <= brd_t)
                        | int(iprevPt_y >= brd_b);
                    if( x != 0 )
                    {
                        if( level == 0 )
                        {
                            status[ptidx] = 0;
                            if( error )
                                error[ptidx] = 0;
                        }
                        continue;
                    }

                    var a:Number = prevPt_x - iprevPt_x;
                    var b:Number = prevPt_y - iprevPt_y;
                    const W_BITS14:int = 14;
                    const W_BITS4:int = 14;
                    const W_BITS1m5:int = W_BITS4 - 5;
                    const W_BITS1m51:int = (1 << ((W_BITS1m5) - 1));
                    const W_BITS14_:int = (1 << W_BITS14);
                    const W_BITS41:int = (1 << ((W_BITS4) - 1));
                    
                    const FLT_SCALE:Number = 1.0/(1 << 20);

                    var iw00:int = ((1.0 - a)*(1.0 - b)*W_BITS14_) + 0.5;
                    var iw01:int = (a*(1.0 - b)*W_BITS14_) + 0.5;
                    var iw10:int = ((1.0 - a)*b*W_BITS14_) + 0.5;
                    var iw11:int = __cint(W_BITS14_ - iw00 - iw01 - iw10);

                    var step:int = lev_w;
                    var A11:Number = 0, A12:Number = 0, A22:Number = 0;
                    var A11i:int = 0, A12i:int = 0, A22i:int = 0;

                    // extract the patch from the first image, compute covariation matrix of derivatives
                    for( y = 0; y < winSize; ++y )
                    {
                        var src:int = __cint( (y + iprevPt_y)*step + iprevPt_x );
                        var dsrc:int = src << 1;//__cint( (y + iprevPt_y)*dstep + iprevPt_x*2 );

                        var Iptr:int = __cint(y*IWinBuf_step);
                        var dIptr:int = __cint(y*derivIWinBuf_step);

                        x = 0;
                        for( ; x < winSize; )
                        {
                            var ival:int = __cint( (img_prev[src])*iw00 + (img_prev[src+1])*iw01 +
                                    (img_prev[src+step])*iw10 + (img_prev[src+step+1])*iw11 );
                            ival = __cint(((ival) + W_BITS1m51) >> (W_BITS1m5));

                            var ixval:int = __cint( deriv[dsrc]*iw00 + deriv[dsrc+2]*iw01 +
                                    deriv[dsrc+dstep]*iw10 + deriv[dsrc+dstep+2]*iw11 );
                            ixval = __cint(((ixval) + W_BITS41) >> (W_BITS4));

                            var iyval:int = __cint( deriv[dsrc+1]*iw00 + deriv[dsrc+3]*iw01 + deriv[dsrc+dstep+1]*iw10 +
                                    deriv[dsrc+dstep+3]*iw11 );
                            iyval = __cint(((iyval) + W_BITS41) >> (W_BITS4));

                            IWinBuf[Iptr] = ival;
                            derivIWinBuf[dIptr] = ixval;
                            __asm(IncLocalInt(Iptr),IncLocalInt(dIptr));
                            derivIWinBuf[dIptr] = iyval;
                            __asm(IncLocalInt(dIptr));

                            A11i = __cint(A11i +ixval*ixval);
                            A12i = __cint(A12i +ixval*iyval);
                            A22i = __cint(A22i +iyval*iyval);

                            dsrc = __cint(dsrc + 2);
                            __asm(IncLocalInt(src),IncLocalInt(x));
                        }
                    }

                    A11 = Number(A11i) * FLT_SCALE;
                    A12 = Number(A12i) * FLT_SCALE;
                    A22 = Number(A22i) * FLT_SCALE;

                    var D:Number = A11*A22 - A12*A12;
                    var minEig:Number = (A22 + A11 - Math.sqrt((A11-A22)*(A11-A22) +
                                        4.0*A12*A12)) / Number(winArea2);

                    if( null != error && (flags & OPTFLOW_GET_MIN_EIGENVALS) != 0 )
                        error[ptidx] = minEig;

                    if( minEig < minEigThreshold || D < FLT_EPSILON )
                    {
                        if( level == 0 )
                            status[ptidx] = 0;
                        continue;
                    }

                    D = 1.0/D;

                    nextPt_x -= halfWin;
                    nextPt_y -= halfWin;
                    prevDelta_x = 0.;
                    prevDelta_y = 0.;

                    for( j = 0; j < maxIterations; ++j )
                    {
                        inextPt_x = (nextPt_x);
                        inextPt_y = (nextPt_y);
                        
                        x = int(inextPt_x <= brd_l) 
                            | int(inextPt_x >= brd_r) 
                            | int(inextPt_y <= brd_t)
                            | int(inextPt_y >= brd_b);
                        if( x != 0 )
                        {
                            if( level == 0 )
                                status[ptidx] = 0;
                            break;
                        }

                        a = nextPt_x - inextPt_x;
                        b = nextPt_y - inextPt_y;
                        iw00 = ((1.0 - a)*(1.0 - b)*W_BITS14_) + 0.5;
                        iw01 = (a*(1.0 - b)*W_BITS14_) + 0.5;
                        iw10 = ((1.0 - a)*b*W_BITS14_) + 0.5;
                        iw11 = __cint(W_BITS14_ - iw00 - iw01 - iw10);
                        var b1:Number = 0, b2:Number = 0;
                        var b1i:int = 0, b2i:int = 0;

                        for( y = 0; y < winSize; )
                        {
                            var Jptr:int = __cint( (y + inextPt_y)*step + inextPt_x );

                            Iptr = __cint(y*IWinBuf_step);
                            dIptr = __cint(y*derivIWinBuf_step);

                            x = 0;
                            for( ; x < winSize; )
                            {
                                var diff:int = __cint( (img_next[Jptr])*iw00 + (img_next[Jptr+1])*iw01 +
                                        (img_next[Jptr+step])*iw10 + (img_next[Jptr+step+1])*iw11 );
                                diff = __cint(((diff) + W_BITS1m51) >> (W_BITS1m5));
                                diff = __cint(diff - IWinBuf[Iptr]);

                                b1i = __cint(b1i + diff*derivIWinBuf[dIptr]);
                                __asm(IncLocalInt(dIptr));
                                b2i = __cint(b2i + diff*derivIWinBuf[dIptr]);
                                __asm(IncLocalInt(dIptr));

                                __asm(IncLocalInt(Iptr));
                                __asm(IncLocalInt(Jptr));
                                __asm(IncLocalInt(x));
                            }
                            __asm(IncLocalInt(y));
                        }

                        b1 = Number(b1i) * FLT_SCALE;
                        b2 = Number(b2i) * FLT_SCALE;

                        var delta_x:Number = ((A12*b2 - A22*b1) * D);
                        var delta_y:Number = ((A12*b1 - A11*b2) * D);

                        nextPt_x += delta_x;
                        nextPt_y += delta_y;
                        nextPts[id0] = nextPt_x + halfWin;
                        nextPts[id1] = nextPt_y + halfWin;

                        if( delta_x*delta_x + delta_y*delta_y <= epsilon )
                            break;

                        if( j > 0 && FastMath.abs(delta_x + prevDelta_x) < 0.01 &&
                                FastMath.abs(delta_y + prevDelta_y) < 0.01 )
                        {
                            nextPts[id0] -= delta_x*0.5;
                            nextPts[id1] -= delta_y*0.5;
                            break;
                        }
                        prevDelta_x = delta_x;
                        prevDelta_y = delta_y;
                    } // iterations loop

                    // calc error
                    if( null != error && level == 0 && status[ptidx] && (flags & OPTFLOW_GET_MIN_EIGENVALS) == 0 )
                    {
                        nextPoint_x = nextPts[id0] - halfWin;
                        nextPoint_y = nextPts[id1] - halfWin;

                        inextPoint_x = (nextPoint_x);
                        inextPoint_y = (nextPoint_y);
                        
                        x = int(inextPoint_x <= brd_l) 
                            | int(inextPoint_x >= brd_r) 
                            | int(inextPoint_y <= brd_t)
                            | int(inextPoint_y >= brd_b);

                        if( x != 0 )
                        {
                            status[ptidx] = 0;
                            continue;
                        }

                        var aa:Number = nextPoint_x - inextPoint_x;
                        var bb:Number = nextPoint_y - inextPoint_y;
                        iw00 = ((1.0 - aa)*(1.0 - bb)*W_BITS14_) + 0.5;
                        iw01 = (aa*(10. - bb)*W_BITS14_) + 0.5;
                        iw10 = ((1.0 - aa)*bb*W_BITS14_) + 0.5;
                        iw11 = __cint(W_BITS14_ - iw00 - iw01 - iw10);
                        //var errval:Number = 0.0;
                        
                        // changed to compute NCC
                        var sum0:int = 0, sum1:int = 0;
                        var sum_sq0:int = 0, sum_sq1:int = 0;
                        var sum_01:int = 0;
                        for( y = 0; y < winSize; )
                        {
                            Jptr = __cint( (y + inextPoint_y)*step + inextPoint_x );
                            Iptr = __cint( y*IWinBuf_step );

                            for( x = 0; x < winSize; )
                            {
                                /*
                                diff = __cint( (img_next[Jptr])*iw00 + (img_next[Jptr+1])*iw01 +
                                        (img_next[Jptr + step] ) * iw10 + (img_next[Jptr + step + 1] ) * iw11 );
                                diff = __cint(((diff) + W_BITS41) >> (W_BITS4));
                                diff = __cint(diff - IWinBuf[Iptr]);
                                errval += IntMath.abs(diff);
                                */
                                ival = __cint( (img_next[Jptr])*iw00 + (img_next[Jptr+1])*iw01 +
                                        (img_next[Jptr + step] ) * iw10 + (img_next[Jptr + step + 1] ) * iw11 );
                                ival = __cint(((ival) + W_BITS1m51) >> (W_BITS1m5));
                                
                                ixval = IWinBuf[Iptr];
                                
                                sum0 = __cint(sum0 + ival);
                                sum1 = __cint(sum1 + ixval);

                                sum_sq0 = __cint(sum_sq0 + (ival*ival));
                                sum_sq1 = __cint(sum_sq1 + (ixval*ixval));
                                sum_01 = __cint(sum_01 + (ival*ixval));
                                
                                __asm(IncLocalInt(Jptr), IncLocalInt(Iptr));
                                __asm(IncLocalInt(x));
                            }
                            __asm(IncLocalInt(y));
                        }
                        var sum0f:Number = Number(sum0) * FLT_SCALE; 
                        var sum1f:Number = Number(sum1) * FLT_SCALE;
                        var sum_sq0f:Number = Number(sum_sq0) * FLT_SCALE; 
                        var sum_sq1f:Number = Number(sum_sq1) * FLT_SCALE;
                        var sum_01f:Number = Number(sum_01) * FLT_SCALE;

                        var num:Number = sum_01f - invWin * sum1f * sum0f;

                        error[ptidx] = (num*num) / ((sum_sq1f - invWin * sum1f * sum1f) * (sum_sq0f - invWin * sum0f * sum0f));
                        //error[ptidx] = errval * 1.0/(32*winArea);
                    }
                } // points loop
            } // levels loop
        }
    }
}
