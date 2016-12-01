classdef MinHeap < handle
    %MINHEAP Min Heap class
    %   copyright (c) 2016 Zorah Lähner (laehner@in.tum.de)
    
    properties
        maxSize_;
        heap_;
        map_;
        size_;
        minEmpty_;
        maxFull_;
    end
    
    methods
        function obj = MinHeap(maxSize)
            obj.maxSize_ = maxSize;
            obj.heap_ = -1 * ones(maxSize,2);
            obj.map_ = -1 * ones(maxSize,1);
            obj.size_ = 0;
            obj.minEmpty_ = 1;
            obj.maxFull_ = 0;
        end
        
        function [obj, b] = decrease(obj, key, newValue)
            if (key < 1 || key > obj.maxSize_)
               b = false;
               return;
            end
            
            index = obj.map_(key);
            % check if bubbling is needed
            if (index ~= -1 && newValue < obj.heap_(index,1))
                obj.heap_(index,1) = newValue;
                obj = obj.bubbleUp(index);
                % check if position actually changed
                if (obj.heap_(index,1) ~= newValue) 
                   b = true; 
                   return;
                end
            end
            
            b = false;
        end
        
        function [obj, b] = increase(obj, key, newValue)
            if(key < 1 || key > obj.maxSize_)
               b = false;
               return;
            end
            
            index = obj.map_(key);
            % check if bubbling is needed
            if (index ~= -1 && newValue > obj.heap_(index,1)) 
               obj.heap_(index,1) = newValue;
               obj = obj.bubbleDown(index);
               % check if position actually changed
               if(obj.heap_(index,1) ~= newValue)
                  b = true;
                  return;
               end
            end
            
            b = false;
        end
        
        function b = isEmpty(obj) 
           b = (obj.size_ == 0); 
        end
        
        function min = peak(obj) 
            if (obj.size_ == 0), min = Inf; else min = obj.heap_(1,1); end
        end
        
        function value = peakKey(obj, key)
            if (key < 1 || key > obj.maxSize_ || obj.map_(key) == -1)
               value = -1;
               return;
            end
            
            value = obj.heap_(obj.map_(key),1);
        end
        
        function [obj, minimum] = pop(obj)
            if (obj.size_ == 0)
               minimum = [-1, -1]; 
               return;
            end
            
            minimum = obj.heap_(1,:);
            
            % bubble down to maintain heap property
            obj.heap_(1,1) = Inf;
            obj = obj.switchEntries(1, obj.maxFull_);
            % clean minimum object
            obj.heap_(obj.maxFull_, :) = [-1, -1];
            obj.map_(minimum(1,2)) = -1;
            
            [obj, newIndex] = obj.bubbleDown(1);
            
            % update minimum/maximum empty/full index if necessary
            if(obj.maxFull_ < obj.minEmpty_)
               obj.minEmpty_ = obj.maxFull_; 
            end
            
            while(obj.maxFull_ >= 1 && obj.heap_(obj.maxFull_, 2) == -1)
               obj.maxFull_ = obj.maxFull_ - 1; 
            end
            
            if(newIndex > obj.maxFull_)
               obj.maxFull_ = newIndex; 
            end
            
            obj.size_ = obj.size_ - 1;
        end
        
        function obj = push(obj, newValue, key)
            if (key < 1 || key > obj.maxSize_ || obj.map_(key) ~= -1 || obj.size_ == obj.maxSize_)
               return; 
            end
            
            obj.heap_(obj.minEmpty_, :) = [newValue, key];
            obj.map_(key) = obj.minEmpty_;
            obj = obj.bubbleUp(obj.minEmpty_);
            obj.size_ = obj.size_ + 1;
            
            % update minimum/maximum empty/full index
            if (obj.minEmpty_ == obj.maxFull_ + 1) 
               obj.maxFull_ = obj.maxFull_ + 1; 
            end
            
            if (obj.size_ == obj.maxSize_)
               obj.minEmpty_ = obj.size_ + 1;
               return;
            end
            
            while(obj.heap_(obj.minEmpty_, 2) > -1)
               obj.minEmpty_ = obj.minEmpty_ + 1; 
            end
        end
        
        function [obj, i] = bubbleDown(obj, i)
            % i invalid
            if (i < 1 || i > obj.maxSize_)
               i = -1;
               return;
            end
            
            if (i <= 0), left = -1; else left = (i * 2); end
            if (i <= 0), right = -1; else right = (i * 2) + 1; end
            
            while ((left <= obj.maxSize_ || right <= obj.maxSize_) && (left ~= -1 || right ~= -1))
               goLeft = false;
               goRight = false;
               
               leftPossible = (left ~= -1 && left <= obj.maxSize_ && obj.heap_(left,2) ~= -1 && obj.heap_(left,1) <= obj.heap_(i,1));
               rightPossible = (right ~= -1 && right <= obj.maxSize_ && obj.heap_(right,2) ~= -1 && obj.heap_(right,1) <= obj.heap_(i,1));
               
               if(leftPossible && rightPossible)
                   if(obj.heap_(right,1) <= obj.heap_(left,1))
                      goRight = true;
                   else
                       goLeft = true;
                   end
               elseif (leftPossible)
                   goLeft = true;
               elseif (rightPossible)
                   goRight = true;
               end
               
               % bubble process
               if(goLeft)
                   % bubble to the left
                   newParent = obj.heap_(left,:);
                   obj.heap_(left,:) = obj.heap_(i,:);
                   obj.heap_(i,:) = newParent;
                   
                   % update map
                   obj.map_(newParent(1,2)) = i;
                   obj.map_(obj.heap_(left,2)) = left;
                   
                   % update indices
                   i = left;
                   if (i<=0), left = -1; else left = (i*2); end
                   if (i<=0), right = -1; else right = (i*2) + 1; end
               elseif(goRight)
                   % bubble to the right
                   newParent = obj.heap_(right,:);
                   obj.heap_(right,:) = obj.heap_(i,:);
                   obj.heap_(i,:) = newParent;
                   
                   % update map
                   obj.map_(newParent(1,2)) = i;
                   obj.map_(obj.heap_(right,2)) = right;
                   
                   % update indices
                   i = right;
                   if (i<=0), left = -1; else left = (i*2); end
                   if (i<=0), right = -1; else right = (i*2) + 1; end
               else
                   % bubbling ends
                   left = -1;
                   right = -1;
               end
                           
            end
            
        end
        
        function [obj, i] = bubbleUp(obj, i)
            % i invalid
            if (i < 1 || i > obj.maxSize_)
                i = -1;
                return;
            end
            
            if (i<=1), parent = -1; else parent = floor(i/2); end
            while (parent ~= -1)
               if(obj.heap_(parent,1) >= obj.heap_(i,1))
                  % continue bubbling
                  newChild = obj.heap_(parent,:);
                  obj.heap_(parent,:) = obj.heap_(i,:);
                  obj.heap_(i,:) = newChild;
               
                  % update map
                  obj.map_(newChild(1,2)) = i;
                  obj.map_(obj.heap_(parent,2)) = parent;
                  
                  % update indices
                  i = parent;
                  if (i <= 1), parent = -1; else parent = floor(i/2); end
               else
                   % bubbling ends :(
                   parent = -1;
               end
            end
        end
        
        function obj = switchEntries(obj, i, j)
            if(i == j || i < 1 || j < 0 || i > obj.maxSize_ || j > obj.maxSize_ || obj.heap_(i, 2) == -1 || obj.heap_(j,2) == -1)
               return; 
            end
            
            % switch heap
            pairi = obj.heap_(i,:);
            obj.heap_(i,:) = obj.heap_(j,:);
            obj.heap_(j,:) = pairi;
            
            % update map
            obj.map_(pairi(1,2)) = j;
            obj.map_(obj.heap_(i,2)) = i;
        end
    end
    
end

