%%%%其他读正确帧方法，勿删
if size(obj.T_lessframe,2)==0 && size(obj.T_moreframe,2)~=0
                mf=getnonXnum(obj.z_moreframe(find(obj.T_moreframe<ti),:),0);
                if length(find(obj.T_moreframe==ti))==0
                    fi=mf+obj.t0 + obj.T * (ti-1) + zi - 1;
                else
                    mz=getnonXnum(obj.z_moreframe(find(obj.T_moreframe==ti),find(obj.z_moreframe(find(obj.T_moreframe==ti),:)<zi)),0);
                    fi=mf+mz+obj.t0 + obj.T * (ti-1) + zi - 1;
                end
            elseif size(obj.T_lessframe,2)~=0 && size(obj.T_moreframe,2)==0
                lf=getnonXnum(obj.z_lessframe(find(obj.T_lessframe<ti),:),0);
                if length(find(obj.T_lessframe==ti))==0
                    fi=-lf+obj.t0 + obj.T * (ti-1) + zi - 1;
                else
                    lz=getnonXnum(obj.z_lessframe(find(obj.T_lessframe==ti),find(obj.z_lessframe(find(obj.T_lessframe==ti),:)<zi)),0);
                    fi=-lf-lz+obj.t0 + obj.T * (ti-1) + zi - 1;
                end
            elseif size(obj.T_lessframe,2)~=0 && size(obj.T_moreframe,2)~=0
                mf=getnonXnum(obj.z_moreframe(find(obj.T_moreframe<ti),:),0);
                lf=getnonXnum(obj.z_lessframe(find(obj.T_lessframe<ti),:),0);
                if length(find(obj.T_moreframe==ti))~=0
                    mz=getnonXnum(obj.z_moreframe(find(obj.T_moreframe==ti),find(obj.z_moreframe(find(obj.T_moreframe==ti),:)<zi)),0);
                else
                    mz=0;
                end
                if length(find(obj.T_lessframe==ti))~=0
                    lz=getnonXnum(obj.z_lessframe(find(obj.T_lessframe==ti),find(obj.z_lessframe(find(obj.T_lessframe==ti),:)<zi)),0); 
                else
                    lz=0;
                end
                fi=mf+mz-lf-lz+obj.t0 + obj.T * (ti-1) + zi - 1; 
            end