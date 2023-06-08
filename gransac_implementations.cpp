#include "gransac_implementations.h"

void findInliers(vector<NewEvent> &events, NewEvent &processed_event, TBranch *branch_eventnum, TBranch *branch_t, TBranch *branch_chg, TBranch *branch_chmb, TBranch *branch_layer, TBranch *branch_tube, TBranch *branch_is_inliner)
{
    vector<NewEvent> processed_events;

    for (int i = 0; i < events.size(); i++)
    {
        if (i % 5000 == 0)
            std::cout << "you are at event " << i << std::endl;

        std::vector<std::shared_ptr<GRANSAC::AbstractParameter>> points;
        vector<vector<float>> tubecoords;

        for (int j = 0; j < events.at(i).t.size(); j++)
        {
            auto xy = getTubeCoords(events.at(i).chamber.at(j), events.at(i).layer.at(j), events.at(i).tube.at(j));

            tubecoords.push_back(xy);

            auto point = std::make_shared<Point2D>(xy.at(0), xy.at(1));

            points.push_back(point);
        }

        // Set up GRANSAC
        GRANSAC::RANSAC<Line2DModel, 2> gransac;

        // Run GRANSAC

        gransac.Initialize(3, 200);

        gransac.Estimate(points);

        auto BestInliers = gransac.GetBestInliers();

        vector<int> goodhitindexes;

        for (int inliers = 0; inliers < BestInliers.size(); inliers++)
        {
            auto point = std::dynamic_pointer_cast<Point2D>(BestInliers.at(inliers));

            float x = point->m_Point2D[0];
            float y = point->m_Point2D[1];

            for (int j = 0; j < events.at(i).t.size(); j++)
            {
                if (tubecoords.at(j)[0] == x && tubecoords.at(j)[1] == y)
                {
                    goodhitindexes.push_back(j);
                }
            }
        }

        processed_event = events.at(i);

        for (int j = 0; j < events.at(i).t.size(); j++)
        {
            if (std::find(goodhitindexes.begin(), goodhitindexes.end(), j) != goodhitindexes.end())
            {
                events.at(i).is_inlier.at(j) = 1;
                processed_event.is_inlier.at(j) = 1;
            }
            else
            {
                events.at(i).is_inlier.at(j) = 0;
                processed_event.is_inlier.at(j) = 0;
            }
        }

        branch_eventnum->Fill();
        branch_t->Fill();
        branch_chg->Fill();
        branch_chmb->Fill();
        branch_layer->Fill();
        branch_tube->Fill();
        branch_is_inliner->Fill();
    }
}