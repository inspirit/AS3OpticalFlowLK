package  
{
    import flash.display.Bitmap;
    import flash.display.BitmapData;
    import flash.display.Graphics;
    import flash.display.Shape;
	import flash.display.Sprite;
    import flash.display.StageScaleMode;
    import flash.events.Event;
    import flash.events.MouseEvent;
    import flash.filters.BlurFilter;
    import flash.geom.Point;
    import flash.geom.Rectangle;
    import flash.media.Camera;
    import flash.media.Video;
    import flash.text.TextField;
    import flash.ui.ContextMenu;
    import flash.ui.ContextMenuItem;
    import flash.utils.getTimer;
    import ru.inspirit.asfeat.oflow.OpticalFlowPyrLK;
	
	/**
     * ...
     * @author Eugene Zatepyakin
     */
    [SWF(width='640',height='520',frameRate='25',backgroundColor='0xFFFFFF')]
    public final class TestPyrLK extends Sprite 
    {
        public const blur4_filter:BlurFilter = new BlurFilter(4, 4, 1);
		public const blur2_filter:BlurFilter = new BlurFilter(2, 2, 2);
        public const ORIGIN:Point = new Point;
        
        protected var myview:Sprite;
        public static var _txt:TextField;
        protected var camBmp:Bitmap;
        protected var shape:Shape;
        protected var gfx:Graphics;
        
        protected var _cam:Camera;
        protected var _video:Video;
        protected var _cambuff:BitmapData;
        protected var _cambuff_rect:Rectangle;
        
        protected var _buff:BitmapData;
        
        public var curr_pts:Vector.<Number>;
        public var prev_pts:Vector.<Number>;
        public var curr_pts_cnt:int = 0;
        public var pyrLK:OpticalFlowPyrLK;
        
        public var corn_t:int = 0;
        public var corn_tt:int = 0;
        public var corn_tn:int = 0;
        
        public function TestPyrLK() 
        {
            if (stage) init();
            else addEventListener(Event.ADDED_TO_STAGE, init);
        }
        
        protected function init(e:Event = null):void
        {
            removeEventListener(Event.ADDED_TO_STAGE, init);
            initStage();            
            
            myview = new Sprite();
            
            // debug test field
            _txt = new TextField();
            _txt.autoSize = 'left';
            _txt.width = 300;
            _txt.x = 5;
            _txt.y = 480;                   
            myview.addChild(_txt);
            
            // web camera initiation
            initCamera(640, 480, 25);
            
            _cambuff = new BitmapData( 640, 480, false, 0x0 );
			_cambuff_rect = _cambuff.rect;			
			_buff = new BitmapData(640, 480, false, 0x00);
            
            camBmp = new Bitmap(_cambuff);
            myview.addChild(camBmp);
            
            shape = new Shape();
            gfx = shape.graphics;
            myview.addChild(shape);
            
            // allocate enough pts buffer
            curr_pts = new Vector.<Number>(5000 * 2, true);
            prev_pts = curr_pts.concat();
            curr_pts_cnt = 0;
            
            // init pyr LK
            pyrLK = new OpticalFlowPyrLK();
            pyrLK.setup(640, 480, 15, 3, true);
            // init first frame before starting loops
            pyrLK.buildNextImagePyramid(_cambuff, true);
            
            addChild(myview);
			addEventListener(Event.ENTER_FRAME, render);
            stage.addEventListener(MouseEvent.CLICK, onMouseClick);
        }
        
        private function onMouseClick(e:MouseEvent):void 
        {
            var mx:int = myview.mouseX;
            var my:int = myview.mouseY;
            var n:int = curr_pts_cnt * 2;
            var fnd:Boolean = false;
            for (var i:int = 0; i < n; )
            {
                var ix:int = curr_pts[i++];
                var iy:int = curr_pts[i++];
                var dx:int = mx - ix;
                var dy:int = my - iy;
                if (dx * dx + dy * dy < 25)
                {
                    fnd = true;
                    break;
                }
            }
            
            if (fnd)
            {
                var di:int = i - 2;
                for (; i < n; )
                {
                    var nx:Number = curr_pts[i++];
                    var ny:Number = curr_pts[i++];
                    curr_pts[di++] = nx;
                    curr_pts[di++] = ny;
                }
                curr_pts_cnt--;
            } else {
                curr_pts[n++] = mx;
                curr_pts[n] = my;
                curr_pts_cnt++;
            }
        }
        
        public function render(e:Event):void
		{	
			_cambuff.draw( _video );
            
            _txt.text = "";
			
            // swap pts buffers
            var tmp:Vector.<Number> = curr_pts;
            curr_pts = prev_pts;
            prev_pts = tmp;
            
            var t:int = getTimer();
            
            // apply some blur
            _buff.applyFilter(_cambuff, _cambuff_rect, ORIGIN, blur2_filter);
            
            // swap pyrLK buffers
            pyrLK.swapImagePyramids();
            // build pyrLK pyramid
            pyrLK.buildNextImagePyramid(_buff, true);
            // calc pyrLK derivatives
            pyrLK.calcPrevImageDerivatives();
            
            // track points
            var pts_status:Vector.<int> = new Vector.<int>(curr_pts_cnt);
            pyrLK.run(prev_pts, curr_pts, curr_pts_cnt, pts_status, null, 20, 0.01, 0.0001);
            
            corn_t += getTimer() - t;
			corn_tn++;
			if(corn_tn == 25)
			{
				corn_tt = corn_t / corn_tn;
				corn_t = 0;
				corn_tn = 0;
			}
            
            // check results
            var n:int = curr_pts_cnt;
            var j:int = 0, k:int = 0;
            for (var i:int = 0; i < n; ++i, j += 2)
            {
                if (pts_status[i])
                {
                    curr_pts[k++] = curr_pts[j];
                    curr_pts[k++] = curr_pts[(j+1)|0];
                }
            }
            curr_pts_cnt = k >> 1;
			
            gfx.clear();
			plotCurrPoints();
            
            _txt.appendText('\nprocessing time: ' + corn_tt + 'ms');
            _txt.appendText('\npoints count: ' + curr_pts_cnt);
		}
        
        protected function plotCurrPoints():void
		{
			var col:uint = 0x00FF00;
			var px:Number, py:Number;
			var n:int = curr_pts_cnt * 2;
			
            gfx.beginFill(col);
			for(var i:int = 0; i < n; )
			{
				px = curr_pts[i++];
				py = curr_pts[i++];
				
                gfx.drawCircle(px, py, 2);
			}
            gfx.endFill();
		}
        
        protected function initCamera(w:int = 640, h:int = 480, fps:int = 25):void
        {
            _cam = Camera.getCamera();
            _cam.setMode( w, h, fps, true );
            
            _video = new Video( _cam.width, _cam.height );
            _video.attachCamera( _cam );
        }
        
        protected function initStage():void
		{
			stage.scaleMode = StageScaleMode.NO_SCALE;

			var myContextMenu:ContextMenu = new ContextMenu();
			myContextMenu.hideBuiltInItems();

			var copyr:ContextMenuItem;
			copyr = new ContextMenuItem("© inspirit.ru", false, false);
			myContextMenu.customItems.push(copyr);

			contextMenu = myContextMenu;
		}
        
    }

}